"""Recompute, from the stored simulation outputs, every simulation figure the documentation quotes.

    python validation/summarize.py            # print the summary
    python validation/summarize.py --write    # and write validation/summary.json and the
                                              # generated block of validation/README.md

Reads ``validation/data`` (``data/README.md`` says what is there and how it was made) and
the package that is installed; runs no quantifier.  It stops when a stored file is not the
one ``data/README.md`` names by SHA-256.  ``tests/test_validation.py`` holds the summary,
and every figure quoted from it, to what this computes.

For each of the 49 genes the package is asked two things about the sequence the
simulation read: what the v2.1.1 gate said (:func:`old_gate`) and what the structural
system says at the simulated design (:func:`structural`).  The second is the call the
harness made before any read existed (``run_simulation.predicted``), so it has to give
what every stored row recorded; when it does not, the stored comparison no longer
describes the code that ships, and :func:`summarize` says which gene and which field.

Every figure is over one panel: the genes with a row in all seven Salmon designs.  The rank
statistics are the harness's own (``run_simulation.ranks``, ``spearman``,
``mannwhitney``: average ranks, two-sided normal approximation with no continuity
correction).
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import os
import re
import statistics
import sys

import numpy as np

from isoform_dominance import identifiability as I
from isoform_dominance import io

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SUMMARY = os.path.join(HERE, "summary.json")
README = os.path.join(HERE, "README.md")
START, END = "<!-- summary:start -->", "<!-- summary:end -->"


def _harness():
    spec = importlib.util.spec_from_file_location(
        "_run_simulation", os.path.join(HERE, "run_simulation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


H = _harness()

#: bias_5p of the seven designs, 3' skew to 5' skew, and the stored file of each
B = (-2.0, -1.0, -0.3, 0.0, 0.3, 1.0, 2.0)
FILES = {"salmon": ("c1_bm20", "c1_bm10", "c1_bm03", "c1_d1", "c1_b03", "c1_b10", "c1_bias5p"),
         "kallisto": ("k1_bm20", "k1_bm10", "k1_bm03", "k1_d1", "k1_b03", "k1_b10", "k1_b20")}
RERUN = "chk_d1"
UNIFORM = B.index(0.0)
FIVE = [i for i, b in enumerate(B) if b > 0]
THREE = [i for i, b in enumerate(B) if b < 0]
SKEWED = THREE + FIVE

#: the design the simulation predicted at: run_simulation's defaults
DESIGN = {"tpm": 10.0, "depth": 30e6, "mean_efflen": 1500.0}

#: the default --tau of 2.3.0, withdrawn: the conditioning bin edge of the table
TAU = 10.0

#: the 1.23-fold ratio the direction was scored at, and the band the CLI states it in
DIRECTION_FOLD = 1.23

#: the rows of issue #2, by (old gate passes, conditioning bin): 1 and 2 pass the gate and are
#: estimable, well and badly conditioned; 3 fails it and is estimable; 4 fails it and is not
ISSUE_ROWS = {(True, "<= 10"): "1", (True, "> 10"): "2", (False, "<= 10"): "3",
              (False, "> 10"): "3", (False, "inf"): "4"}


# --------------------------------------------------------------------------- #
# the stored files
# --------------------------------------------------------------------------- #
def sha256_table(path):
    """``{file: SHA-256}`` from the tables of a README, outside its ``<details>``.

    A row is one whose first cell is a file name in backticks and whose last is a
    SHA-256 in backticks.
    """
    text = open(path, encoding="utf-8").read().split("<details>")[0]
    return {m.group(1): m.group(2) for m in re.finditer(
        r"^\| `([^`]+)`[^\n]*\| `([0-9a-f]{64})` \|$", text, re.M)}


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_files(data=DATA):
    """The stored files whose bytes are not the ones ``data/README.md`` names, as
    ``[(file, why)]``; empty when every one is."""
    bad = []
    for name, want in sorted(sha256_table(os.path.join(data, "README.md")).items()):
        path = os.path.join(data, name)
        if not os.path.exists(path):
            bad.append((name, "missing"))
        elif file_sha256(path) != want:
            bad.append((name, "SHA-256 %s, not %s" % (file_sha256(path), want)))
    return bad


def load(data=DATA):
    """The inputs bundle and the stored rows: ``(genes, rows, class_sizes)``.

    ``genes`` is the bundle's ``{gene: {"config", "sequences"}}``; ``rows`` maps a
    quantifier to seven ``{gene: row}``, one per design of :data:`B`, and ``"rerun"`` to
    ``chk_d1``'s.
    """
    with gzip.open(os.path.join(data, "sweep116.json.gz"), "rt", encoding="utf-8") as f:
        genes = json.load(f)["genes"]
    rows = {}
    for quant, names in FILES.items():
        rows[quant] = []
        for b, name in zip(B, names, strict=True):
            with open(os.path.join(data, "sim", name + ".json"), encoding="utf-8") as f:
                doc = json.load(f)
            # the file is the design it is named for, quantified by the tool it is filed under
            assert doc["design"]["bias_5p"] == b, name
            assert doc["design"].get("quantifier", "salmon") == quant, name
            rows[quant].append({r["gene"]: r for r in doc["genes"]})
    with open(os.path.join(data, "sim", RERUN + ".json"), encoding="utf-8") as f:
        rows["rerun"] = {r["gene"]: r for r in json.load(f)["genes"]}
    with open(os.path.join(data, "sim", "class_sizes.json"), encoding="utf-8") as f:
        sizes = json.load(f)
    return genes, rows, sizes


# --------------------------------------------------------------------------- #
# what the package says about each gene
# --------------------------------------------------------------------------- #
def _split(entry):
    """``(config, configured sequences, the rest of the gene)``, ids without version,
    as ``run_simulation.predicted`` splits them."""
    cfg = entry["config"]
    seqs = {t.split(".")[0]: s for t, s in entry["sequences"].items()}
    need = {t.split(".")[0] for v in cfg["groups"].values() for t in v}
    return (cfg, {t: seqs[t] for t in seqs if t in need},
            {t: s for t, s in seqs.items() if t not in need})


def old_gate(entry):
    """The v2.1.1 gate: does each class of the primary comparison have a strand-aware
    31-mer that no other configured class has?  The rest of the gene is not consulted."""
    cfg, conf, _ = _split(entry)
    res = I.analyze(cfg, sequences=conf, canonical=False, background_gene_transcripts=False)
    n = {g: res["groups"][g]["n_unique_kmers"] for g in res["primary_comparison"]}
    return {"passes": all(v > 0 for v in n.values()), "n_unique_kmers": n}


def structural(entry):
    """The structural verdict at the simulated design, as the harness asked for it, and
    the whole report, which the survey figures read."""
    cfg, conf, rest = _split(entry)
    res = I.analyze(cfg, sequences=conf, background_sequences=rest,
                    background_gene_transcripts=False, **DESIGN)
    c, pc = res["contrast"], res["primary_comparison"]
    out = {"verdict": res["verdict"], "estimable": c["estimable"],
           "conditioning_factor": c["conditioning_factor"],
           "min_resolvable_log2fc": c["min_resolvable_log2fc"],
           "beyond_linear": c["beyond_linear"], "log2_efflen_ratio": c["log2_efflen_ratio"],
           "primary_comparison": list(pc),
           "transcript_lengths": {g: [len(conf[t.split(".")[0]]) for t in cfg["groups"][g]]
                                  for g in pc}}
    # the other estimand, which decides nothing the harness recorded: the fields above are
    # the fragment estimand's, as every stored prediction is
    ab = res["estimands"]["abundance"]
    out["abundance"] = {"estimable": ab["contrast"]["estimable"],
                        "min_resolvable_log2fc": ab["contrast"]["min_resolvable_log2fc"]}
    return out, res


def predicted_se(s):
    """``min_resolvable_log2fc / 1.96``, as the harness records it (``inf`` when there is
    no finite figure)."""
    fc = s["min_resolvable_log2fc"]
    return fc / 1.96 if math.isfinite(fc) else float("inf")


def _same_float(a, b):
    if math.isinf(a) or math.isinf(b):
        return math.isinf(a) and math.isinf(b)
    return math.isclose(a, b, rel_tol=1e-9)


def consistency(structure, rows):
    """Every stored row whose prediction is not what :func:`structural` gives now, as
    ``[(gene, file, field, stored, now)]``, and how many rows were compared."""
    out, n = [], 0
    sources = [(name, d) for quant, names in FILES.items()
               for name, d in zip(names, rows[quant], strict=True)] + [(RERUN, rows["rerun"])]
    for name, d in sources:
        for gene, r in sorted(d.items()):
            s, n = structure[gene], n + 1
            now = {"verdict": s["verdict"], "estimable": s["estimable"],
                   "conditioning": s["conditioning_factor"],
                   "min_resolvable_log2fc": s["min_resolvable_log2fc"],
                   "beyond_linear": s["beyond_linear"], "predicted_se": predicted_se(s)}
            for field, v in now.items():
                same = (_same_float(float(r[field]), float(v)) if isinstance(v, float)
                        else r[field] == v)
                if not same:
                    out.append((gene, name, field, r[field], v))
    return out, n


# --------------------------------------------------------------------------- #
# the figures
# --------------------------------------------------------------------------- #
def _sign(x):
    return (x > 0) - (x < 0)


def _counting_floor(size):
    """SD of log2(nA/nB) from Poisson counting alone, at the drawn class sizes."""
    return math.sqrt(1.0 / size["nA"] + 1.0 / size["nB"]) / math.log(2.0)


def _median(values):
    values = list(values)
    return statistics.median(values) if values else None


def _direction(designs, genes, ratio):
    """Per design, how many genes' bias took the 5' direction (``sign(bias) =
    -sign(ratio)``: the shorter class inflated) and how many the 3' (``+sign(ratio)``)."""
    return {"five_prime_rule": [sum(_sign(d[g]["bias"]) == -_sign(ratio[g]) for g in genes)
                                for d in designs],
            "three_prime_rule": [sum(_sign(d[g]["bias"]) == _sign(ratio[g]) for g in genes)
                                 for d in designs]}


def _efflen(lengths):
    return [I.effective_length(n) for n in lengths]


#: the class-length summaries scored against bias (``class_efflen_ratio``'s docstring):
#: each gives log2 of class A's figure over class B's, positive when class A is longer
LENGTH_SUMMARIES = {
    "plain mean": lambda a, b: math.log2((sum(a) / len(a)) / (sum(b) / len(b))),
    "harmonic mean": lambda a, b: math.log2((len(a) / sum(1 / x for x in a))
                                            / (len(b) / sum(1 / x for x in b))),
    "class total": lambda a, b: math.log2(sum(a) / sum(b)),
    "minimum efflen": lambda a, b: math.log2(min(a) / min(b)),
    # S = sum of 1/efflen is large for a short class, so log2(S_B / S_A)
    "sum of 1/efflen": lambda a, b: math.log2(sum(1 / x for x in b) / sum(1 / x for x in a)),
}


def _poisson_factor(entry, theta_of):
    """``sqrt(c'(A'WA)^+ c)`` of the class contrast with ``W = diag(1 / A theta)``, the
    windows at k = 31, and ``theta`` of each transcript ``theta_of(its sequence)``."""
    cfg, conf, rest = _split(entry)
    groups = {g: [t.split(".")[0] for t in v] for g, v in cfg["groups"].items()}
    pc = cfg["primary_comparison"]
    ids = [t for v in groups.values() for t in v] + sorted(rest)
    seqs = dict(conf, **rest)
    A, _ = I.compatibility_matrix({t: I.kmer_track(seqs[t], I.DEFAULT_K) for t in ids}, ids)
    c = np.array([1.0 if t in groups[pc[0]] else -1.0 if t in groups[pc[1]] else 0.0
                  for t in ids])
    cov = I.gls_covariance(A, np.array([float(theta_of(seqs[t])) for t in ids]))
    return math.sqrt(float(c @ cov @ c))


def poisson_unit_factor(entry):
    """The Poisson-weighted counterpart of the conditioning factor that
    ``gls_covariance``'s docstring compares it with, at ``theta = 1`` per transcript."""
    return _poisson_factor(entry, lambda seq: 1.0)


def design_fragments(seq):
    """The expected fragments of a transcript at the simulated design (:data:`DESIGN`):
    ``identifiability.expected_informative_reads`` with every fragment counted."""
    return I.expected_informative_reads(1.0, DESIGN["tpm"], len(seq), depth=DESIGN["depth"],
                                        mean_efflen=DESIGN["mean_efflen"])


def poisson_design_factor(entry):
    """The same factor with the design's expected fragment counts as ``theta``
    (:func:`design_fragments`), which differ from one transcript to another with its
    length: the other figure ``gls_covariance``'s docstring gives."""
    return _poisson_factor(entry, design_fragments)


def analyse(genes):
    """What the package says about each gene: ``{gene: {"old_gate", "structural",
    "coherence", "counting_min_resolvable_log2fc", "poisson_unit_factor",
    "poisson_design_factor"}}``.  The slow part of :func:`summarize`, and a function of the
    sequence alone."""
    out = {}
    for gene in sorted(genes):
        s, res = structural(genes[gene])
        out[gene] = {
            "old_gate": old_gate(genes[gene]), "structural": s,
            "coherence": {c: res["groups"][c]["coherence"]["min_jaccard"]
                          for c in res["primary_comparison"]},
            "counting_min_resolvable_log2fc":
                res["unique_read_counting_se"]["min_resolvable_log2fc"],
            "poisson_unit_factor": poisson_unit_factor(genes[gene]) if s["estimable"] else None,
            "poisson_design_factor":
                poisson_design_factor(genes[gene]) if s["estimable"] else None}
    return out


def summarize(data=DATA, analysed=None):
    """The whole summary, as the dict ``summary.json`` holds (non-finite numbers as
    ``inf``; :func:`write` turns them into null).  ``analysed`` is :func:`analyse`'s answer
    for this data's sequence, when the caller has it."""
    bad = check_files(data)
    if bad:
        raise SystemExit("validation/data is not what data/README.md names:\n"
                         + "\n".join("  %s: %s" % b for b in bad))
    genes, rows, sizes = load(data)
    salmon, kallisto = rows["salmon"], rows["kallisto"]
    analysed = analysed if analysed is not None else analyse(genes)
    structure = {g: analysed[g]["structural"] for g in genes}
    per = {g: {"old_gate": analysed[g]["old_gate"], "structural": structure[g],
               "rows": {q: sum(g in d for d in rows[q]) for q in FILES}} for g in sorted(genes)}
    mismatches, n_rows = consistency(structure, rows)

    # ---- the panel --------------------------------------------------------- #
    panel = [g for g in sorted(genes) if all(g in d for d in salmon)]
    # kallisto's, said beside it: S7 reads kallisto's rows of the Salmon panel
    kallisto_panel = [g for g in sorted(genes) if all(g in d for d in kallisto)]
    for g in per:
        per[g]["in_panel"] = g in panel
    left_out = {g: {"salmon_designs_with_row": [b for b, d in zip(B, salmon, strict=True)
                                                 if g in d],
                    "kallisto_designs_with_row": [b for b, d in zip(B, kallisto, strict=True)
                                                   if g in d]}
                for g in sorted(genes) if g not in panel}
    est = [g for g in panel if structure[g]["estimable"]
           and math.isfinite(predicted_se(structure[g]))]
    ni = [g for g in panel if not structure[g]["estimable"]]
    ratio = {g: structure[g]["log2_efflen_ratio"] for g in panel}
    floor = {g: _counting_floor(sizes[g]) for g in panel}

    def excess(d, g):
        return d[g]["rmse"] / floor[g]

    def worst(g):
        return max(abs(d[g]["bias"]) for d in salmon)

    # ---- S1: spread, not accuracy ----------------------------------------- #
    s1 = {"genes": est, "n": len(est),
          "abs_bias_over_predicted_se": [sum(abs(d[g]["bias"]) > d[g]["predicted_se"]
                                             for g in est) for d in salmon],
          "realised_sd_under_predicted_se": [sum(d[g]["realised_sd"] < d[g]["predicted_se"]
                                                 for g in est) for d in salmon],
          "median_rmse_over_predicted_se": [statistics.median(d[g]["rmse"] / d[g]["predicted_se"]
                                                              for g in est) for d in salmon],
          "NTRK3": {"in_S1": "NTRK3" in est,
                    "realised_sd": [min(d["NTRK3"]["realised_sd"] for d in salmon),
                                    max(d["NTRK3"]["realised_sd"] for d in salmon)],
                    "abs_bias": [min(abs(d["NTRK3"]["bias"]) for d in salmon),
                                 max(abs(d["NTRK3"]["bias"]) for d in salmon)]}}

    # ---- S2, S7: the direction, Salmon and kallisto ------------------------ #
    band = [g for g in panel if abs(ratio[g]) >= I.DIRECTION_MIN_ABS_LOG2_EFFLEN_RATIO]
    fold = [g for g in panel if abs(ratio[g]) > math.log2(DIRECTION_FOLD)]
    s2 = dict(_direction(salmon, band, ratio), genes=band, n=len(band),
              not_identifiable=[g for g in band if g in ni],
              same_genes_as_over_1_23_fold=band == fold)
    s7 = dict(_direction(kallisto, band, ratio), n=len(band))

    # ---- S3: LEPR ------------------------------------------------------------ #
    want = {i: -_sign(ratio["LEPR"]) if B[i] > 0 else _sign(ratio["LEPR"]) for i in SKEWED}
    s3 = {"log2_efflen_ratio": ratio["LEPR"], "skewed_designs": len(SKEWED),
          "against_rule": sum(_sign(salmon[i]["LEPR"]["bias"]) == -want[i] for i in SKEWED),
          "bias": [d["LEPR"]["bias"] for d in salmon]}

    # ---- S4: which class-length summary ----------------------------------- #
    s4 = {"n": len(panel), "summaries": {}}
    for name, f in LENGTH_SUMMARIES.items():
        r = {g: f(*(_efflen(structure[g]["transcript_lengths"][c])
                    for c in structure[g]["primary_comparison"])) for g in panel}
        big = [g for g in panel if abs(r[g]) > math.log2(DIRECTION_FOLD)]
        s4["summaries"][name] = {
            "spearman": [H.spearman([r[g] for g in panel], [d[g]["bias"] for g in panel])
                         for d in salmon],
            "five_prime_sign_right": sum(_sign(salmon[B.index(2.0)][g]["bias"]) == -_sign(r[g])
                                         for g in big),
            "n_over_1_23_fold": len(big)}

    # ---- S5: does not_identifiable flag a larger bias ---------------------- #
    s5 = {"n_not_identifiable": len(ni), "n_estimable": len(est), "not_identifiable": ni,
          "mannwhitney_p": [H.mannwhitney([abs(d[g]["bias"]) for g in ni],
                                          [abs(d[g]["bias"]) for g in est])[1]
                            for d in salmon]}

    # ---- S6: does conditioning, or the predicted SE, rank the error -------- #
    s6 = {"n": len(est),
          "conditioning_vs_excess": [H.spearman(
              [structure[g]["conditioning_factor"] for g in est], [excess(d, g) for g in est])
              for d in salmon],
          "predicted_se_vs_rmse": [H.spearman(
              [d[g]["predicted_se"] for g in est], [d[g]["rmse"] for g in est])
              for d in salmon]}

    # ---- S8: the design's confound ----------------------------------------- #
    def ratio_of(f):
        a, b = _median(f(g) for g in ni), _median(f(g) for g in est)
        return {"not_identifiable": a, "estimable": b, "ratio": a / b}
    s8 = {"median_gene_fragments": ratio_of(lambda g: sizes[g]["total"]),
          "median_smaller_class_fragments": ratio_of(
              lambda g: min(sizes[g]["nA"], sizes[g]["nB"])),
          "counting_floor_mannwhitney_p": H.mannwhitney([floor[g] for g in ni],
                                                        [floor[g] for g in est])[1]}

    # ---- the four-row table of issue #2 ------------------------------------- #
    table = []
    for gate in (True, False):
        for bin_ in ("<= 10", "> 10", "inf"):
            gs = [g for g in panel if per[g]["old_gate"]["passes"] == gate
                  and (not structure[g]["estimable"] if bin_ == "inf" else
                       structure[g]["estimable"]
                       and (structure[g]["conditioning_factor"] <= TAU) == (bin_ == "<= 10"))]
            table.append({"old_gate": "pass" if gate else "fail",
                          "structure": "not estimable" if bin_ == "inf" else "estimable",
                          "conditioning": bin_, "n": len(gs),
                          "issue_row": ISSUE_ROWS.get((gate, bin_)),
                          "uniform_rmse_median": _median(salmon[UNIFORM][g]["rmse"] for g in gs),
                          "uniform_excess_median": _median(excess(salmon[UNIFORM], g)
                                                           for g in gs),
                          "worst_abs_bias_median": _median(worst(g) for g in gs),
                          "genes": gs if len(gs) <= 10 else None})

    def row(gate, bin_):
        return [g for g in panel if per[g]["old_gate"]["passes"] == gate
                and structure[g]["estimable"]
                and (structure[g]["conditioning_factor"] <= TAU) == (bin_ == "<= 10")]
    good, high = row(True, "<= 10"), row(True, "> 10")
    comparison = {
        "uniform_rmse": H.mannwhitney([salmon[UNIFORM][g]["rmse"] for g in high],
                                      [salmon[UNIFORM][g]["rmse"] for g in good])[1],
        "uniform_excess": H.mannwhitney([excess(salmon[UNIFORM], g) for g in high],
                                        [excess(salmon[UNIFORM], g) for g in good])[1],
        "worst_abs_bias": H.mannwhitney([worst(g) for g in high], [worst(g) for g in good])[1]}
    notes = {"ESR1": {"worst_abs_bias": worst("ESR1"),
                      "replicates_with_a_ratio": [d["ESR1"]["n"] for d in salmon]},
             "NTRK3_uniform_rmse": salmon[UNIFORM]["NTRK3"]["rmse"],
             "estimable_uniform_rmse_median": _median(salmon[UNIFORM][g]["rmse"] for g in est)}

    # ---- the 49-gene survey the code's docstrings quote --------------------- #
    have = sorted(g for g in genes if structure[g]["estimable"])
    cond = sorted(structure[g]["conditioning_factor"] for g in have)
    unit = {g: analysed[g]["poisson_unit_factor"] for g in have}
    at_design = {g: analysed[g]["poisson_design_factor"] for g in have}
    coherence = sorted((j, g, c) for g in sorted(genes)
                       for c, j in analysed[g]["coherence"].items() if j is not None)
    lengths = [len(s) for g in genes.values() for s in g["sequences"].values()]
    survey = {"n_genes": len(genes), "n_estimable": len(have),
              "n_transcripts": len(lengths), "shortest_transcript": min(lengths),
              "median_conditioning": statistics.median(cond),
              "conditioning_range": cond[-1] / cond[0],
              "fraction_over_tau": sum(x > TAU for x in cond) / len(cond),
              "poisson_unit_median": statistics.median(unit.values()),
              "poisson_unit_range": max(unit.values()) / min(unit.values()),
              "poisson_unit_spearman": H.spearman([structure[g]["conditioning_factor"]
                                                   for g in have], [unit[g] for g in have]),
              "poisson_design_median": statistics.median(at_design.values()),
              "poisson_design_range": max(at_design.values()) / min(at_design.values()),
              "poisson_design_spearman": H.spearman(
                  [structure[g]["conditioning_factor"] for g in have],
                  [at_design[g] for g in have]),
              "classes_under_0_05_min_jaccard": [{"gene": g, "class": c, "min_jaccard": j}
                                                 for j, g, c in coherence if j < 0.05],
              # over the primary-comparison classes with two or more members, the others
              # having no pair: class_coherence's docstring
              "n_classes_with_pairs": len(coherence),
              "median_min_jaccard": statistics.median(j for j, _, _ in coherence),
              "TPI1_counting_min_resolvable_log2fc":
                  analysed["TPI1"]["counting_min_resolvable_log2fc"],
              # the contrast as the abundance estimand: estimable for how many, and where the
              # two estimands part, which genes, each way
              "n_estimable_abundance": sum(structure[g]["abundance"]["estimable"]
                                           for g in genes),
              "estimable_as_abundance_only": sorted(
                  g for g in genes if structure[g]["abundance"]["estimable"]
                  and not structure[g]["estimable"]),
              "estimable_as_fragment_only": sorted(
                  g for g in genes if structure[g]["estimable"]
                  and not structure[g]["abundance"]["estimable"])}

    return {"format": "isoform-dominance/validation-summary/1",
            "designs": {"bias_5p": list(B), "files": {q: list(v) for q, v in FILES.items()},
                        "rerun": RERUN, "predicted_at": DESIGN},
            "genes": per,
            "consistency": {"rows_compared": n_rows,
                            "mismatches": [list(m) for m in mismatches]},
            "panel": {"rule": "a row in all seven Salmon designs", "genes": panel,
                      "n": len(panel), "left_out": left_out,
                      "kallisto_genes": kallisto_panel, "kallisto_n": len(kallisto_panel)},
            "S1": s1, "S2": s2, "S3": s3, "S4": s4, "S5": s5, "S6": s6, "S7": s7, "S8": s8,
            "table": table, "table_comparison": comparison, "table_notes": notes,
            "survey": survey}


# --------------------------------------------------------------------------- #
# the generated block of validation/README.md
# --------------------------------------------------------------------------- #
def _f(x, places):
    return "" if x is None else "%.*f" % (places, x)


def _span(values, places, dash="–"):
    lo, hi = min(values), max(values)
    return _f(lo, places) if _f(lo, places) == _f(hi, places) else (
        _f(lo, places) + dash + _f(hi, places))


def render(s):
    """The Markdown between the summary markers of validation/README.md."""
    sal = s["S1"]
    bs = "| design | " + " | ".join("%+.1f" % b if b else "0" for b in B) + " |"
    sep = "|---|" + "---|" * len(B)
    out = []
    w = out.append
    w("Generated by `python validation/summarize.py --write` from `validation/data`; "
      "`tests/test_validation.py` fails when it no longer agrees.")
    w("")
    p = s["panel"]
    if p["kallisto_genes"] == p["genes"]:
        kallisto = "kallisto has a row for the same %d" % p["kallisto_n"]
    else:
        kallisto = ("kallisto has a row in all seven of its designs for %d genes: not %s, and "
                    "also %s" % (p["kallisto_n"],
                                 ", ".join(sorted(set(p["genes"]) - set(p["kallisto_genes"])))
                                 or "none",
                                 ", ".join(sorted(set(p["kallisto_genes"]) - set(p["genes"])))
                                 or "none"))
    w("Panel: %d genes, the ones with a row in all seven Salmon designs (%s). Left out: %s."
      % (p["n"], kallisto, "; ".join(
          "%s (Salmon rows at bias_5p %s)" % (g, ", ".join("%+.1f" % b if b else "0"
                                                          for b in v["salmon_designs_with_row"]))
          for g, v in sorted(s["panel"]["left_out"].items()))))
    w("")
    w("Every stored row's prediction (verdict, estimable, conditioning, "
      "min_resolvable_log2fc, beyond_linear, predicted_se; %d rows of the 15 files) is what "
      "the installed package gives for that gene now: %d differ."
      % (s["consistency"]["rows_compared"], len(s["consistency"]["mismatches"])))
    w("")
    w("### The four rows of issue #2")
    w("")
    w("Old gate: the v2.1.1 rule, a strand-aware 31-mer of each class that no other "
      "configured class has. Structure: the class contrast estimable at the simulated design. "
      "Conditioning: the contrast's conditioning factor against %g, the default `--tau` "
      "withdrawn in 2.3.0. RMSE and excess at uniform coverage (excess = RMSE over the "
      "Poisson counting floor of the drawn class sizes); worst |bias|: per gene, the largest "
      "of the seven Salmon designs. Medians over the genes of a row." % TAU)
    w("")
    w("| issue row | old gate | structure | conditioning | n | RMSE | excess | worst \\|bias\\| "
      "| genes |")
    w("|---|---|---|---|---|---|---|---|---|")
    for r in s["table"]:
        w("| %s | %s | %s | %s | %d | %s | %s | %s | %s |" % (
            r["issue_row"] or "(none)", r["old_gate"], r["structure"],
            "∞" if r["conditioning"] == "inf" else r["conditioning"].replace("<=", "≤"),
            r["n"], _f(r["uniform_rmse_median"], 3), _f(r["uniform_excess_median"], 2),
            _f(r["worst_abs_bias_median"], 2), ", ".join(r["genes"] or [])))
    c, n = s["table_comparison"], s["table_notes"]
    w("")
    w("Pass, estimable, conditioning > %g against ≤ %g (Mann-Whitney p): uniform RMSE %s, "
      "uniform excess %s, worst |bias| %s. ESR1: worst |bias| %s; replicates with a defined "
      "ratio per design %s. NTRK3's uniform RMSE: %s; the estimable genes' median: %s."
      % (TAU, TAU, _f(c["uniform_rmse"], 2), _f(c["uniform_excess"], 2),
         _f(c["worst_abs_bias"], 2), _f(n["ESR1"]["worst_abs_bias"], 2),
         "/".join(str(x) for x in n["ESR1"]["replicates_with_a_ratio"]),
         _f(n["NTRK3_uniform_rmse"], 3), _f(n["estimable_uniform_rmse_median"], 3)))
    w("")
    w("### S1. Spread, not accuracy (Salmon; the %d estimable genes with a finite predicted "
      "SE)" % sal["n"])
    w("")
    w(bs)
    w(sep)
    w("| \\|bias\\| > predicted SE | " + " | ".join(
        str(x) for x in sal["abs_bias_over_predicted_se"]) + " |")
    w("| replicate SD < predicted SE | " + " | ".join(
        str(x) for x in sal["realised_sd_under_predicted_se"]) + " |")
    w("| median RMSE / predicted SE | " + " | ".join(
        _f(x, 2) for x in sal["median_rmse_over_predicted_se"]) + " |")
    w("")
    t = sal["NTRK3"]
    w("predicted SE = `min_resolvable_log2fc` / 1.96. NTRK3 (not estimable, so not among "
      "the %d): replicate SD %s, |bias| %s log2 across the seven designs."
      % (sal["n"], _span(t["realised_sd"], 4), _span(t["abs_bias"], 4)))
    for key, label in (("S2", "S2. Direction, Salmon"), ("S7", "S7. Direction, kallisto "
                                                                "(the same reads)")):
        d = s[key]
        w("")
        w("### %s (the %d genes with |log2 effective-length ratio| ≥ %g)"
          % (label, d["n"], I.DIRECTION_MIN_ABS_LOG2_EFFLEN_RATIO))
        w("")
        w(bs)
        w(sep)
        w("| 5' rule: sign(bias) = −sign(ratio) | " + " | ".join(
            str(x) for x in d["five_prime_rule"]) + " |")
        w("| 3' rule: sign(bias) = +sign(ratio) | " + " | ".join(
            str(x) for x in d["three_prime_rule"]) + " |")
    w("")
    w("The %d include %d that are not estimable (%s), and are the same genes as those whose "
      "ratio exceeds %g-fold: %s. kallisto has no bias model and leans to the longer class at "
      "uniform coverage too, so its 3' figures are no evidence of the skew; only the 5' ones "
      "turn toward it."
      % (s["S2"]["n"], len(s["S2"]["not_identifiable"]), ", ".join(s["S2"]["not_identifiable"]),
         DIRECTION_FOLD, "yes" if s["S2"]["same_genes_as_over_1_23_fold"] else "no"))
    l3 = s["S3"]
    w("")
    w("### S3. LEPR")
    w("")
    w("log2 effective-length ratio %s; its bias ran against the rule in %d of the %d skewed "
      "designs (bias by design: %s)." % (_f(l3["log2_efflen_ratio"], 3), l3["against_rule"],
                                         l3["skewed_designs"],
                                         ", ".join(_f(x, 3) for x in l3["bias"])))
    w("")
    w("### S4. Class-length summaries (Salmon, the %d panel genes)" % s["S4"]["n"])
    w("")
    w("| summary | Spearman with bias, b = +2.0 | b = −2.0 | b = 0 | 5' sign right at b = +2.0, "
      "ratio > %g-fold |" % DIRECTION_FOLD)
    w("|---|---|---|---|---|")
    for name, v in s["S4"]["summaries"].items():
        sp = v["spearman"]
        w("| %s | %+.3f | %+.3f | %+.3f | %d/%d |" % (name, sp[B.index(2.0)], sp[B.index(-2.0)],
                                                     sp[UNIFORM], v["five_prime_sign_right"],
                                                     v["n_over_1_23_fold"]))
    w("")
    w("\"sum of 1/efflen\" is log2(S_B/S_A), positive when class A is longer, like the others.")
    w("")
    w("### S5. |bias|, the %d not estimable against the %d estimable (Mann-Whitney p)"
      % (s["S5"]["n_not_identifiable"], s["S5"]["n_estimable"]))
    w("")
    w(bs)
    w(sep)
    w("| p | " + " | ".join(_f(x, 2) for x in s["S5"]["mannwhitney_p"]) + " |")
    w("")
    w("Not estimable: %s. Read S8 first." % ", ".join(s["S5"]["not_identifiable"]))
    w("")
    w("### S6. What ranks the error (Spearman, the %d estimable genes)" % s["S6"]["n"])
    w("")
    w(bs)
    w(sep)
    w("| conditioning vs excess | " + " | ".join(
        "%+.3f" % x for x in s["S6"]["conditioning_vs_excess"]) + " |")
    w("| predicted SE vs RMSE | " + " | ".join(
        "%+.3f" % x for x in s["S6"]["predicted_se_vs_rmse"]) + " |")
    e = s["S8"]
    w("")
    w("### S8. The design's confound")
    w("")
    w("Every transcript has TPM 10, so a gene's fragments grow with its transcript count, "
      "which is also what makes a class contrast inestimable. Not estimable against "
      "estimable, medians: gene fragments %s against %s (%.1fx), fragments of the smaller "
      "class %s against %s (%.1fx); the counting floor differs at Mann-Whitney p = %s."
      % ("{:,.0f}".format(e["median_gene_fragments"]["not_identifiable"]),
         "{:,.0f}".format(e["median_gene_fragments"]["estimable"]),
         e["median_gene_fragments"]["ratio"],
         "{:,.0f}".format(e["median_smaller_class_fragments"]["not_identifiable"]),
         "{:,.0f}".format(e["median_smaller_class_fragments"]["estimable"]),
         e["median_smaller_class_fragments"]["ratio"], _f(e["counting_floor_mannwhitney_p"], 4)))
    v = s["survey"]
    w("")
    w("### The structural survey of the same %d genes" % v["n_genes"])
    w("")
    w("Quoted by the package's docstrings, the `--min-log2fc` help and the harness. %d "
      "transcripts, the shortest %d nt. The contrast is estimable "
      "for %d; their conditioning factor: median %s, range %.0fx, above %g for %.0f%%. "
      "The Poisson-weighted factor at theta = 1 per transcript: median %s, range %.0fx, "
      "Spearman %s with the conditioning factor; with the design's expected fragment "
      "counts as theta: median %s, range %sx, Spearman %s. Classes whose least similar pair of "
      "transcripts shares under 5%% of windows: %s; over the %d classes of the primary "
      "comparisons with two or more transcripts, the median share of the least similar "
      "pair: %s. TPI1's unique-read counting min |log2FC|: %s. As the abundance estimand "
      "(length-corrected class abundance, `--estimand abundance`) the contrast is estimable "
      "for %d: for %d genes as abundance only (%s), and as fragment totals only for %s."
      % (v["n_transcripts"], v["shortest_transcript"],
         v["n_estimable"], _f(v["median_conditioning"], 1), v["conditioning_range"], TAU,
         100 * v["fraction_over_tau"], _f(v["poisson_unit_median"], 1),
         v["poisson_unit_range"], _f(v["poisson_unit_spearman"], 2),
         _f(v["poisson_design_median"], 2), _f(v["poisson_design_range"], 1),
         _f(v["poisson_design_spearman"], 3),
         ", ".join("%s %s %s" % (x["gene"], x["class"], _f(x["min_jaccard"], 3))
                   for x in v["classes_under_0_05_min_jaccard"]),
         v["n_classes_with_pairs"], _f(v["median_min_jaccard"], 2),
         _f(v["TPI1_counting_min_resolvable_log2fc"], 2),
         v["n_estimable_abundance"], len(v["estimable_as_abundance_only"]),
         ", ".join(v["estimable_as_abundance_only"]) or "none",
         ", ".join(v["estimable_as_fragment_only"]) or "none"))
    return "\n".join(out)


def as_json(summary):
    """``summary.json``'s text: keys sorted, a number that is not finite as null."""
    return json.dumps(io.finite_json(summary), indent=1, sort_keys=True,
                      allow_nan=False, ensure_ascii=False) + "\n"


def with_block(text, block):
    """``text`` with the part between the summary markers replaced by ``block``."""
    head, rest = text.split(START, 1)
    _, tail = rest.split(END, 1)
    return "%s%s\n%s\n%s%s" % (head, START, block, END, tail)


def write(summary, path=SUMMARY, readme=README):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(as_json(summary))
    with open(readme, encoding="utf-8") as f:
        text = f.read()
    with open(readme, "w", encoding="utf-8", newline="\n") as f:
        f.write(with_block(text, render(summary)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true",
                    help="write validation/summary.json and the block of validation/README.md")
    a = ap.parse_args(argv)
    s = summarize()
    print(render(s))
    if s["consistency"]["mismatches"]:
        print("\nThe stored predictions are not what the package gives now:", file=sys.stderr)
        for m in s["consistency"]["mismatches"][:20]:
            print("  %s %s %s: stored %r, now %r" % tuple(m), file=sys.stderr)
    if a.write:
        write(s)
        print("\nwrote %s and the block of %s" % (SUMMARY, README), file=sys.stderr)
    return 1 if s["consistency"]["mismatches"] else 0


if __name__ == "__main__":
    sys.exit(main())
