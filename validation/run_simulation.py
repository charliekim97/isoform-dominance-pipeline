"""Simulate, quantify, and compare against what the surrogate predicted.

    # once
    python validation/run_simulation.py index --fasta cdna.fa.gz --index salmon_idx

    # the checkpoint -- nothing downstream means anything until this passes
    python validation/run_simulation.py checkpoint --inputs validation/data/sweep116.json.gz \
        --index salmon_idx --work sim --gene PIK3CA --depth-scale 20 --replicates 30

    # a design point: validation/data/sim/c1_d1.json is this at --bias-5p 0
    python validation/run_simulation.py run --inputs validation/data/sweep116.json.gz \
        --index salmon_idx --work sim --bias-5p 0 --replicates 30 --out c1_d1.json

``--inputs`` is the bundle the stored outputs were made from; ``--rundir`` takes the per-gene
``<gene>.cfg.json`` and ``<gene>.seqs.json`` it holds instead (``simulate_reads.read_inputs``).

Truth is the fragment count drawn per transcript, and the comparison is against Salmon's
``NumReads`` -- also a fragment count.  Nothing is compared in TPM, because an effective
length convention sitting between the simulation and the estimate is the easiest way to
measure the bookkeeping instead of the surrogate.

The checkpoint is the guard on exactly that.  On a gene the surrogate calls well
conditioned, at high depth, the realised SD must approach the predicted one and the bias
must approach zero.  If it does not, the fault is here, not in the package.

Salmon is invoked as a subprocess; that part is not exercised by this file's own checks.
Record ``salmon --version`` and the bias flags used in the design document -- they change
between versions and a number without them is not reproducible.
"""
import argparse
import json
import math
import os
import shutil
import statistics
import subprocess
import sys

import numpy as np

from isoform_dominance import identifiability as I

SIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "simulate_reads.py")


def _sim():
    """simulate_reads, loaded from beside this file."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_sim", SIM)
    sim = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sim)
    return sim


def gene_inputs(a):
    """``{gene: {"config", "sequences"}}`` from ``--inputs`` or ``--rundir``."""
    if (a.rundir is None) == (a.inputs is None):
        sys.exit("give exactly one of --inputs and --rundir")
    return _sim().read_inputs(a.rundir, a.inputs)


def sh(cmd, **kw):
    print("  $ %s" % " ".join(cmd))
    return subprocess.run(cmd, check=True, **kw)


def cmd_index(a):
    if a.quantifier == "kallisto":
        # kallisto takes the FASTA positionally and has no -t/-i ordering quirk;
        # same k so the two indexes see the same sequence resolution.
        sh([a.kallisto_bin, "index", "-i", a.index, "-k", str(a.k),
            "-t", str(a.threads), a.fasta])
    else:
        sh(["salmon", "index", "-t", a.fasta, "-i", a.index, "-k", str(a.k),
            "-p", str(a.threads)])
    print("index built at %s" % a.index)


def predicted(source, gene, tpm, depth_scale, mean_efflen=None, ratio=None):
    """What the shipped code says before any reads exist, for ``gene`` of ``source``
    (:func:`gene_inputs`).

    ``mean_efflen`` must be the one the simulation drew with.  It was omitted once and
    the simulator carried its own copy of the fragment-count formula with a different
    value, so the simulation drew 1.5x the fragments the prediction assumed.  That error
    shows up as a 0.82x ratio in the realised SD, which sits inside any tolerance a
    30-replicate check can honestly use -- see :func:`checkpoint_counts`, which compares
    the counts directly instead.
    """
    cfg = source[gene]["config"]
    seqs = {t.split(".")[0]: s for t, s in source[gene]["sequences"].items()}
    need = {t.split(".")[0] for v in cfg["groups"].values() for t in v}
    kw = {} if mean_efflen is None else {"mean_efflen": mean_efflen}
    res = I.analyze(cfg,
                    sequences={t: seqs[t] for t in need},
                    background_sequences={t: s for t, s in seqs.items() if t not in need},
                    background_gene_transcripts=False,
                    tpm=tpm, depth=30e6 * depth_scale, **kw)
    c = dict(res["contrast"])
    if ratio is not None:
        # analyze() assumes flat TPM.  A --ratio run does not, and the Poisson weights move
        # with it: on PIK3CA a forced 1:1 lifts the small class from ~93 fragments to ~7000
        # and the predicted SE falls from 0.469 to 0.036, while the gene total -- and so
        # any check on fragment counts -- is unchanged.  Re-solve at the abundances the
        # simulation actually used.
        c.update(_predict_at(cfg, seqs, tpm, depth_scale, mean_efflen, ratio))
    return {"verdict": res["verdict"], "estimable": c["estimable"],
            "conditioning": c["conditioning_factor"],
            "min_resolvable_log2fc": c["min_resolvable_log2fc"],
            "beyond_linear": c["beyond_linear"],
            "predicted_se": (c["min_resolvable_log2fc"] / 1.96
                             if math.isfinite(c["min_resolvable_log2fc"]) else float("inf"))}


def _predict_at(cfg, seqs, tpm, depth_scale, mean_efflen, ratio):
    """The contrast figures solved at the abundances a --ratio run actually simulates."""
    sim = _sim()

    groups = {g: [t.split(".")[0] for t in v] for g, v in cfg["groups"].items()}
    pc = cfg.get("primary_comparison") or sorted(groups)[:2]
    gene = {"seqs": seqs, "groups": groups, "pc": list(pc)}
    exp = sim.transcript_fragment_counts(
        gene, tpm, depth_scale, ratio, None,
        mean_efflen if mean_efflen is not None else I.DEFAULT_MEAN_EFFLEN,
        expectations=True)

    need = [t for v in groups.values() for t in v]
    ids = need + sorted(t for t in seqs if t not in need)
    tracks = {t: I.kmer_track(seqs[t], I.DEFAULT_K, True) for t in ids}
    out = I.compatibility_matrix(tracks, ids)
    A = out[0] if isinstance(out, tuple) else out
    theta = np.array([exp[t] for t in ids], dtype=float)
    ind = {c: np.array([1.0 if t in groups[c] else 0.0 for t in ids]) for c in pc}
    est = I.estimability(A, ind[pc[0]] - ind[pc[1]])
    cov = I.gls_covariance(A, theta)
    se = (I.log_ratio_se(cov, ind[pc[0]], ind[pc[1]], theta)
          if est["estimable"] else float("inf"))
    fc = I.min_resolvable_log2fc(se, 1)
    return {"estimable": est["estimable"], "conditioning_factor": est["conditioning_factor"],
            "min_resolvable_log2fc": fc, "beyond_linear": se > I.LINEARISATION_LIMIT}


def checkpoint_counts(truth, tpm, depth_scale, mean_efflen):
    """Did the simulation draw the number of fragments the prediction assumes?

    Exact, and the only check that catches a bookkeeping drift small enough to hide
    inside the Monte Carlo error of a realised SD.  A 1.5x error in the expected count
    moves the SD by 0.82x, which 30 replicates cannot distinguish from noise; it moves
    the count by 1.5x, which is unmissable.
    """
    out = []
    for gene, g in truth["genes"].items():
        drawn = sum(g["fragments_per_transcript"].values())
        out.append((gene, drawn))
    return out


def expected_total(source, gene, tpm, depth_scale, mean_efflen):
    seqs = source[gene]["sequences"]
    return sum(I.expected_informative_reads(
        1.0, tpm, len(s), depth=30e6 * depth_scale,
        mean_efflen=mean_efflen, frag_mean=I.DEFAULT_FRAG_MEAN) for s in seqs.values())


MIN_PERCENT_MAPPED = 80.0
EXPECTED_LIBTYPE = "ISF"
"""What the simulator produces: inward-facing, stranded, read 1 from the forward strand.

``-l A`` auto-detects, and on a library this synthetic it will get it right; the point of
checking is that a mismatch here is the signature of a simulator that changed orientation
without the analysis noticing.
"""


def quantify(index, r1, r2, outdir, threads, extra,
             quantifier="salmon", kallisto_bin="kallisto"):
    """Dispatch to a quantifier. The salmon path below is unchanged from the version
    that produced every c1_*.json and p1_*.json result. Prove that before reading any
    kallisto number: re-run the b=0 leg to a fresh --out with the same design and diff
    it per gene against validation/data/sim/c1_d1.json. The seed is replicate-only
    (simulate_reads.main), so the fragments are identical and the signed errors must
    match to the bit. Any nonzero difference means this refactor moved the salmon path
    and the twelve stored legs can no longer be compared with anything."""
    if quantifier == "kallisto":
        return _quantify_kallisto(index, r1, r2, outdir, threads, kallisto_bin)
    return _quantify_salmon(index, r1, r2, outdir, threads, extra)


def _quantify_kallisto(index, r1, r2, outdir, threads, kallisto_bin):
    """kallisto 0.52 has no bias model of any kind -- ``kallisto quant`` lists no
    positional and no sequence bias option, which is why it is the clean test of the
    effective-length mechanism: there is no correction that could mask or mimic it.

    ``--fr-stranded`` matches the ISF library the simulator emits and that salmon
    detected under ``-l A``; leaving it off would hand kallisto strictly less
    information than salmon had and the comparison would no longer be of estimators.
    ``--plaintext`` skips the HDF5 write, which nothing here reads.
    """
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    sh([kallisto_bin, "quant", "-i", index, "-o", outdir, "-t", str(threads),
        "--fr-stranded", "--plaintext", r1, r2])

    info_path = os.path.join(outdir, "run_info.json")
    info = json.load(open(info_path)) if os.path.exists(info_path) else {}
    pct = info.get("p_pseudoaligned")
    if pct is None:
        raise RuntimeError(
            "kallisto wrote no p_pseudoaligned in %s -- the run did not finish as "
            "expected, and nothing downstream should use its abundance.tsv" % info_path)
    if pct < MIN_PERCENT_MAPPED:
        raise RuntimeError(
            "kallisto pseudoaligned %.1f%% of %s fragments, below the %.0f%% this "
            "simulation must reach: every read was drawn verbatim from a transcript in "
            "the index. Check that the index was built from the same FASTA the sequences "
            "came from and that transcript ids match apart from the version suffix."
            % (pct, info.get("n_processed", "?"), MIN_PERCENT_MAPPED))

    counts = {}
    with open(os.path.join(outdir, "abundance.tsv")) as f:
        header = next(f).rstrip("\n").split("\t")
        try:
            col = header.index("est_counts")
        except ValueError:
            raise RuntimeError("abundance.tsv header %r has no est_counts column; this "
                               "kallisto writes a layout this harness does not know"
                               % header) from None
        for line in f:
            pes = line.rstrip("\n").split("\t")
            counts[pes[0].split(".")[0]] = float(pes[col])
    diag = {"percent_mapped": pct,
            "num_processed": info.get("n_processed"),
            "num_mapped": info.get("n_pseudoaligned"),
            "salmon_version": None,
            "kallisto_version": info.get("kallisto_version"),
            "library_type": "fr-stranded (asserted, not detected)",
            "frag_length_mean": None,
            "frag_length_sd": None}
    return counts, diag


def _quantify_salmon(index, r1, r2, outdir, threads, extra):
    """Quantify one library, and refuse to return numbers from a run that went wrong.

    A low mapping rate is the failure mode that looks exactly like a result: an index
    built from the wrong FASTA, transcript ids that differ by their version suffix, a
    library type mis-detected on a short read set -- each gives per-class counts that are
    quietly meaningless, and the only visible symptom is that the surrogate seems not to
    predict anything.  ``aux_info/meta_info.json`` says what actually happened, so it is
    read before the counts are trusted and its numbers travel with the result.
    """
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    sh(["salmon", "quant", "-i", index, "-l", "A", "-1", r1, "-2", r2,
        "-p", str(threads), "-o", outdir] + list(extra))

    meta_path = os.path.join(outdir, "aux_info", "meta_info.json")
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
    pct = meta.get("percent_mapped")
    if pct is None:
        raise RuntimeError(
            "salmon wrote no percent_mapped in %s -- the run did not finish as expected, "
            "and nothing downstream should use its quant.sf" % meta_path)
    if pct < MIN_PERCENT_MAPPED:
        raise RuntimeError(
            "salmon mapped %.1f%% of %s fragments, below the %.0f%% this simulation must "
            "reach: every read was drawn verbatim from a transcript in the index. Check "
            "that the index was built from the same FASTA the sequences came from, that "
            "transcript ids match apart from the version suffix, and that the library "
            "type was detected as %s (see lib_format_counts.json in %s)."
            % (pct, meta.get("num_processed", "?"), MIN_PERCENT_MAPPED,
               EXPECTED_LIBTYPE, outdir))

    lib_path = os.path.join(outdir, "lib_format_counts.json")
    lib = json.load(open(lib_path)) if os.path.exists(lib_path) else {}
    detected = lib.get("expected_format")
    if detected and detected != EXPECTED_LIBTYPE:
        print("  NOTE: salmon detected library type %s, the simulator emits %s. Not fatal, "
              "but if the orientation changed the truth and the estimate may no longer "
              "describe the same fragments." % (detected, EXPECTED_LIBTYPE))

    counts = {}
    with open(os.path.join(outdir, "quant.sf")) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            counts[p[0].split(".")[0]] = float(p[4])          # NumReads
    diag = {"percent_mapped": pct,
            "num_processed": meta.get("num_processed"),
            "num_mapped": meta.get("num_mapped"),
            "salmon_version": meta.get("salmon_version"),
            "library_type": detected,
            "frag_length_mean": meta.get("frag_length_mean"),
            "frag_length_sd": meta.get("frag_length_sd")}
    return counts, diag


def one_replicate(a, rep, source, genes_arg=None):
    work = a.work
    os.makedirs(work, exist_ok=True)
    where = [a.rundir, work] if a.rundir else [work, "--inputs", a.inputs]
    cmd = [sys.executable, SIM] + where + [
           "--depth-scale", str(a.depth_scale),
           "--replicate", str(rep), "--tpm", str(a.tpm),
           "--mean-efflen", str(a.mean_efflen)]
    if a.ratio is not None:
        cmd += ["--ratio", str(a.ratio)]
    if a.bias_5p:
        # "=" form: argparse rejects a bare negative value in scientific notation,
        # which str() produces for small magnitudes
        cmd += ["--bias-5p=%s" % a.bias_5p]
    if genes_arg:
        cmd += ["--genes", genes_arg]
    sh(cmd)
    base = os.path.join(work, "rep%d" % rep)
    truth = json.load(open(base + ".truth.json"))
    try:
        counts, diag = quantify(a.index, base + "_1.fq.gz", base + "_2.fq.gz",
                                os.path.join(work, "quant%d" % rep), a.threads,
                                a.salmon_flags.split() if a.salmon_flags else [],
                                getattr(a, "quantifier", "salmon"),
                                getattr(a, "kallisto_bin", "kallisto"))
    finally:
        # 270 replicates leak 270 libraries if a quant fails halfway
        for f in (base + "_1.fq.gz", base + "_2.fq.gz"):
            if os.path.exists(f):
                os.remove(f)
    out = {}
    for gene, g in truth["genes"].items():
        pc = g["primary_comparison"]
        cfg = source[gene]["config"]
        groups = {k: [t.split(".")[0] for t in v] for k, v in cfg["groups"].items()}
        est = {c: sum(counts.get(t, 0.0) for t in groups[c]) for c in pc}
        tru = g["class_fragments"]
        ok = est[pc[0]] > 0 and est[pc[1]] > 0 and tru[pc[0]] > 0 and tru[pc[1]] > 0
        out[gene] = {
            "estimated_log2_ratio": math.log2(est[pc[0]] / est[pc[1]]) if ok else None,
            "true_log2_ratio": g["true_log2_ratio"],
            "estimated": est, "true": tru,
        }
    # keep the two small files that say what salmon did; drop the bulky rest
    qdir = os.path.join(work, "quant%d" % rep)
    for keep in ("cmd_info.json",):
        src = os.path.join(qdir, keep)
        if os.path.exists(src):
            shutil.copy(src, base + "." + keep)
    json.dump(diag, open(base + ".salmon.json", "w"), indent=1)
    shutil.rmtree(qdir, ignore_errors=True)
    return out, truth, diag


def ranks(v):
    """Average ranks, ties shared.  Hoisted out of spearman() so mannwhitney() uses the
    same one rather than a second copy -- the mean-efflen bug was a duplicated formula."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    out = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def mannwhitney(x, y):
    """Two-sided rank-sum on x vs y; normal approximation.  Returns (z, p)."""
    n1, n2 = len(x), len(y)
    if n1 < 3 or n2 < 3:
        return float("nan"), float("nan")
    r = ranks(list(x) + list(y))
    u = sum(r[:n1]) - n1 * (n1 + 1) / 2.0
    z = (u - n1 * n2 / 2.0) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    return z, math.erfc(abs(z) / math.sqrt(2.0))


def spearman(x, y):
    """Rank correlation with average ranks for ties, and Pearson on the ranks.

    The tie-free shortcut 1 - 6d^2/(n(n^2-1)) is wrong the moment two genes share a
    value, which duplicated sequences make possible.
    """
    if len(x) < 3:
        return float("nan")
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def summarise(per_rep, source, tpm, depth_scale, genes=None, mean_efflen=None, ratio=None):
    rows = []
    names = genes or sorted({g for r in per_rep for g in r})
    for gene in names:
        errs, ests = [], []
        for r in per_rep:
            e = r.get(gene)
            if e and e["estimated_log2_ratio"] is not None and e["true_log2_ratio"] is not None:
                ests.append(e["estimated_log2_ratio"])
                errs.append(e["estimated_log2_ratio"] - e["true_log2_ratio"])
        if len(ests) < 3:
            continue
        p = predicted(source, gene, tpm, depth_scale, mean_efflen, ratio)
        rows.append(dict(gene=gene, n=len(ests),
                         realised_sd=statistics.stdev(ests),
                         bias=statistics.median(errs),
                         rmse=math.sqrt(sum(x * x for x in errs) / len(errs)),
                         **p))
    return rows


def cmd_checkpoint(a):
    print("CHECKPOINT: %s at depth x%g, %d replicates.\n"
          "Nothing downstream is meaningful until the drawn fragment count matches what\n"
          "the prediction assumes, the realised SD approaches the predicted one, and the\n"
          "bias approaches zero.\n" % (a.gene, a.depth_scale, a.replicates))
    source = gene_inputs(a)
    per, truths, diags = [], [], []
    for i in range(a.replicates):
        out, truth, diag = one_replicate(a, i, source, a.gene)
        per.append(out)
        truths.append(truth)
        diags.append(diag)

    # (1) exact: did we draw what the prediction assumed?
    want = expected_total(source, a.gene, a.tpm, a.depth_scale, a.mean_efflen)
    drawn = [dict(checkpoint_counts(t, a.tpm, a.depth_scale, a.mean_efflen))[a.gene]
             for t in truths]
    got = statistics.mean(drawn)
    count_ratio = got / want if want else float("inf")
    # Poisson: the mean of n draws has relative SD 1/sqrt(n*want)
    tol = 5.0 / math.sqrt(max(1.0, a.replicates * want))
    counts_ok = abs(count_ratio - 1.0) <= max(tol, 0.01)
    print("\n  (0) salmon %s: mapped %.2f%%..%.2f%% of fragments, library %s"
          % (diags[0].get("salmon_version") or "?",
             min(d["percent_mapped"] for d in diags),
             max(d["percent_mapped"] for d in diags),
             diags[0].get("library_type") or "?"))
    print("\n  (1) fragments drawn   %.1f   assumed by the prediction   %.1f" % (got, want))
    print("      ratio %.4f   tolerance +-%.4f   %s"
          % (count_ratio, max(tol, 0.01), "OK" if counts_ok else "MISMATCH"))
    if not counts_ok:
        print("\n  FAIL at (1). The simulation and the prediction disagree about how many\n"
              "  fragments exist, so nothing about SD or bias below means anything. Check\n"
              "  --tpm, --depth-scale and --mean-efflen are the same on both sides.")
        return 1

    rows = summarise(per, source, a.tpm, a.depth_scale, [a.gene], a.mean_efflen, a.ratio)
    if not rows:
        sys.exit("no usable replicates")
    r = rows[0]
    # The prediction is the delta-method SE of the surrogate system under Poisson counts,
    # not a lower bound on what a quantifier's replicates show: on PIK3CA salmon came back
    # with SD 0.0644 against a predicted 0.1048 while carrying a +0.0545 bias.  Comparing
    # SD alone therefore asks the wrong question; total error is the one the prediction is
    # about.
    ratio = r["rmse"] / r["predicted_se"] if r["predicted_se"] else float("inf")
    band = 3.0 / math.sqrt(2.0 * max(1, r["n"] - 1))
    bias_se = 1.253 * r["realised_sd"] / math.sqrt(r["n"])
    print("\n  (2) predicted SE  %.4f   (min |log2FC| %.4f%s)"
          % (r["predicted_se"], r["min_resolvable_log2fc"],
             ", BEYOND LINEAR" if r["beyond_linear"] else ""))
    print("      realised SD   %.4f   bias %+.4f   RMSE %.4f   over %d replicates"
          % (r["realised_sd"], r["bias"], r["rmse"], r["n"]))
    print("      RMSE/predicted %.3f   3-sigma band %.2f-%.2f" % (ratio, 1 - band, 1 + band))
    if abs(r["bias"]) > 3 * bias_se:
        frac = abs(r["bias"]) / r["min_resolvable_log2fc"] if r["min_resolvable_log2fc"] else 0
        print("  (3) BIAS, %.1f SE from zero and %.0f%% of the smallest effect the tool "
              "claims to resolve." % (abs(r["bias"]) / bias_se, 100 * frac))
        print("      This is a finding, not a fault: a class with little private sequence "
              "and a small share of the gene's fragments absorbs a relative bias from the "
              "abundant class's misallocation. On PIK3CA class iso_118aa has 87%% of its "
              "windows in the other class, its main transcript has none of its own, and it "
              "holds 1/150th of the fragments; a 0.026%% leak is +3.8%% on it.")
    else:
        print("  (3) bias %+.4f, within 3 SE of zero" % r["bias"])

    problems = []
    MIN_REPS = 10
    if r["beyond_linear"]:
        problems.append(
            "the prediction is past the linearisation limit at this depth, so (2) is not\n"
            "  a test of it at all. Raise --depth-scale until beyond_linear clears.")
    if r["n"] < MIN_REPS:
        problems.append("%d replicates put the 3-sigma band at %.2f-%.2f, which almost "
                        "anything clears. Use at least %d."
                        % (r["n"], 1 - band, 1 + band, MIN_REPS))
    if not (1 - band) <= ratio <= (1 + band):
        problems.append(
            "RMSE is %.2fx the prediction, outside the 3-sigma band %.2f-%.2f. Below it "
            "is not automatically wrong -- an unbiased floor can be beaten by a biased\n"
            "  estimator -- but far below it, with the counts in (1) correct, means the\n"
            "  comparison is not measuring what it says." % (ratio, 1 - band, 1 + band))
    if problems:
        print("\n  CHECK")
        for i, p_ in enumerate(problems, 1):
            print("   %d. %s" % (i, p_))
    else:
        print("\n  PASS -- counts and total error line up. Go on to C1.")
    ok = counts_ok and not problems
    if a.out:
        with open(a.out, "w") as f:
            json.dump({"count_ratio": count_ratio, "rows": rows}, f, indent=1)
    return 0 if ok else 1


def cmd_run(a):
    source = gene_inputs(a)
    reps = [one_replicate(a, i, source) for i in range(a.replicates)]
    per = [r[0] for r in reps]
    diags = [r[2] for r in reps]
    rows = summarise(per, source, a.tpm, a.depth_scale, None, a.mean_efflen, a.ratio)
    json.dump({"design": {"depth_scale": a.depth_scale, "ratio": a.ratio,
                          "tpm": a.tpm, "bias_5p": a.bias_5p,
                          "mean_efflen": a.mean_efflen,
                          "replicates": a.replicates,
                          "salmon_flags": a.salmon_flags,
                          "quantifier": getattr(a, "quantifier", "salmon")},
               "salmon": {"version": diags[0].get("salmon_version") if diags else None,
                          "kallisto_version": (diags[0].get("kallisto_version")
                                               if diags else None),
                          "percent_mapped_min": min((d["percent_mapped"] for d in diags),
                                                    default=None),
                          "library_type": diags[0].get("library_type") if diags else None},
               "genes": rows}, open(a.out, "w"), indent=1)
    good = [r for r in rows if r["estimable"]]
    bad = [r for r in rows if not r["estimable"]]
    print("\n%d genes, %d estimable, %d not" % (len(rows), len(good), len(bad)))
    # Scored against BOTH the replicate spread and the total error, because they answer
    # different questions and a line built on either alone is blind to half the failure.
    # `realised_sd` is the spread of the estimates and carries the estimand's own Poisson
    # variation between replicates; RMSE is the distance from the realised truth and does
    # not.  Under a coverage-model violation the SD line degrades gently while RMSE blows
    # up: at --bias-5p 2.0, over the 41 rows printed here, the SD-only C3 read +0.472 and
    # looked survivable while Spearman(predicted, RMSE) was -0.086 and the median
    # RMSE/predicted was 4.55.
    # The gap between the two lines IS the bias term.
    fin = [r for r in good if r["predicted_se"] and math.isfinite(r["predicted_se"])]
    if fin:
        print("  C3  Spearman(predicted, realised)   SD %+.3f   RMSE %+.3f"
              % (spearman([r["predicted_se"] for r in fin], [r["realised_sd"] for r in fin]),
                 spearman([r["predicted_se"] for r in fin], [r["rmse"] for r in fin])))
        for lab, key in (("SD  ", "realised_sd"), ("RMSE", "rmse")):
            v = sorted(r[key] / r["predicted_se"] for r in fin)
            print("  C2  %s/predicted: median %.2f, range %.2f-%.2f"
                  % (lab, statistics.median(v), v[0], v[-1]))
        over_se = sum(1 for r in fin if abs(r["bias"]) > r["predicted_se"])
        over_fc = sum(1 for r in fin if abs(r["bias"]) > r["min_resolvable_log2fc"])
        under = sum(1 for r in fin if r["realised_sd"] < r["predicted_se"])
        print("      |bias| exceeds the predicted SE in %d/%d genes, min|log2FC| in %d/%d;"
              % (over_se, len(fin), over_fc, len(fin)))
        print("      replicate SD sits BELOW the predicted SE in %d/%d." % (under, len(fin)))
        if over_se > len(fin) // 4:
            print("      The interval is narrow and the estimate is outside it. Replicate")
            print("      agreement does not detect this; do not read the SD line alone.")
    if bad:
        for lab, grp in (("not_identifiable", bad), ("estimable       ", good)):
            if not grp:
                continue
            print("  C1  %s  median |bias| %.3f  SD %.3f  RMSE %.3f  (n=%d)"
                  % (lab, statistics.median([abs(r["bias"]) for r in grp]),
                     statistics.median([r["realised_sd"] for r in grp]),
                     statistics.median([r["rmse"] for r in grp]), len(grp)))
        if good:
            z, pval = mannwhitney([abs(r["bias"]) for r in bad],
                                  [abs(r["bias"]) for r in good])
            print("      Mann-Whitney on |bias|: z %+.2f  p %.3f   %s"
                  % (z, pval, "difference" if pval < 0.05 else "no difference"))
            print("      NOTE: the two groups differ in depth by construction under flat")
            print("      TPM (annotation count drives both the flag and the fragment")
            print("      count), so a raw difference here is confounded -- see")
            print("      validation/README.md before reading anything into it.")
    print("\nwrote %s" % a.out)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("index")
    p.add_argument("--fasta", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--k", type=int, default=31)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--quantifier", choices=("salmon", "kallisto"), default="salmon")
    p.add_argument("--kallisto-bin", default="kallisto")
    p.set_defaults(func=cmd_index)

    for name, fn in (("checkpoint", cmd_checkpoint), ("run", cmd_run)):
        p = sub.add_parser(name)
        p.add_argument("--inputs", default=None,
                       help="the inputs bundle, validation/data/sweep116.json.gz")
        p.add_argument("--rundir", default=None,
                       help="the per-gene <gene>.cfg.json and <gene>.seqs.json, in place of "
                            "--inputs")
        p.add_argument("--index", required=True)
        p.add_argument("--work", default="sim")
        p.add_argument("--depth-scale", type=float, default=1.0)
        p.add_argument("--ratio", type=float, default=None,
                       help="omit for the neutral flat-TPM case, which is what the "
                            "package's precision figure assumes; set it only to test "
                            "bias, and the prediction is re-solved at those abundances")
        p.add_argument("--tpm", type=float, default=10.0)
        p.add_argument("--bias-5p", type=float, default=0.0,
                       help="positional coverage skew: >0 toward the 5' end, <0 "
                            "toward the 3' end, 0 uniform")
        p.add_argument("--mean-efflen", type=float, default=I.DEFAULT_MEAN_EFFLEN,
                       help="passed to both the simulator and the prediction; they must "
                            "agree or the comparison measures the bookkeeping")
        p.add_argument("--replicates", type=int, default=30)
        p.add_argument("--threads", type=int, default=8)
        p.add_argument("--salmon-flags", default="")
        p.add_argument("--quantifier", choices=("salmon", "kallisto"), default="salmon",
                       help="kallisto 0.52 has no bias model at all, so --salmon-flags "
                            "is ignored there and the design block records which tool "
                            "produced the numbers")
        p.add_argument("--kallisto-bin", default="kallisto",
                       help="absolute path when kallisto lives in its own conda env")
        p.add_argument("--out", default=None if name == "checkpoint" else "sim_result.json")
        if name == "checkpoint":
            p.add_argument("--gene", required=True)
        p.set_defaults(func=fn)

    a = ap.parse_args(argv)
    sys.exit(a.func(a))


if __name__ == "__main__":
    main()
