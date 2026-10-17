"""Simulate paired-end fragments from real transcript sequences, with the truth kept.

    python simulate_reads.py RUNDIR OUTDIR [--depth-scale 1.0] [--ratio 1.0]
                             [--replicate 0] [--tpm 10] [--genes A,B,...]

Writes OUTDIR/rep<N>_1.fq.gz, rep<N>_2.fq.gz and OUTDIR/rep<N>.truth.json.  The truth is
the number of fragments actually drawn from each transcript -- never a quantifier's own
output, and counted in fragments rather than TPM so that no effective-length convention
sits between the simulation and the comparison.

**Deliberately not the surrogate's model.**  The compatibility matrix is built from whole
transcript window sets; this draws a fragment length, places it somewhere along the
transcript, and emits 100 bp from each end with substitution error.  If a prediction made
from the first still holds against reads made by the second, the two systems were not
told the same thing.  Positional bias is available (`--bias-5p`, negative for the 3' direction) and off by
default, so its effect can be isolated rather than confounded with everything else.

The class ratio is applied to the two classes of each gene's `primary_comparison`, with
the gene's total fragment count held fixed, so variance and bias can be separated: at
1:1 a class total estimate is symmetric, and at 8:1 a quantifier that spreads mass onto
the wrong class shows it as bias rather than only as spread.

**Known design confound, deliberately left in place (2026-09-13).**  TPM is flat *per
transcript*, so a gene's total fragment count is proportional to how many transcripts it
has.  Annotation depth is also what produces a class-crossing null direction, so the
`not_identifiable` group is 4-8x deeper than the estimable one by construction (smaller
class 2611 vs 343 fragments; counting floor differs at Mann-Whitney p = 0.0013).  Every
group comparison between the two is confounded by that, and validation/README.md says so
where it matters.  A future sweep should hold the *gene* total fixed across genes instead.

It is not fixed here because it cannot be fixed on this side alone: `run_simulation.py`
derives its prediction from the package at the same tpm/depth, so rescaling the gene
total in the simulator without re-solving the prediction at the new abundances reproduces
exactly the mean-efflen failure -- a silent constant offset that no checkpoint tolerance
catches.  Doing it properly means a `_predict_at`-style re-solve, which is a change to a
simulator whose validation is closed.  Documented rather than half-done.
"""
import argparse
import glob
import gzip
import json
import math
import os
import sys

import numpy as np

from isoform_dominance import identifiability as I

READ_LEN = 100
FRAG_MEAN, FRAG_SD = I.DEFAULT_FRAG_MEAN, I.DEFAULT_FRAG_SD
ERROR_RATE = 0.002
QUAL = "I"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
BASES = np.array(list("ACGT"))


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def load_genes(rundir, only=None):
    """Every gene the survey has both a config and sequences for."""
    out = {}
    for cfgp in sorted(glob.glob(os.path.join(rundir, "*.cfg.json"))):
        gene = os.path.basename(cfgp)[:-len(".cfg.json")]
        if only and gene not in only:
            continue
        seqp = os.path.join(rundir, "%s.seqs.json" % gene)
        if not os.path.exists(seqp):
            continue
        cfg = json.load(open(cfgp))
        seqs = {t.split(".")[0]: s for t, s in json.load(open(seqp)).items()}
        groups = {g: [t.split(".")[0] for t in v] for g, v in cfg["groups"].items()}
        pc = cfg.get("primary_comparison") or sorted(groups)[:2]
        out[gene] = {"groups": groups, "pc": list(pc), "seqs": seqs}
    return out


def transcript_fragment_counts(gene, tpm, depth_scale, ratio, rng,
                               mean_efflen=I.DEFAULT_MEAN_EFFLEN, expectations=False):
    """Expected fragments per transcript, then a Poisson draw around it.

    The expectation comes from the package's own
    :func:`identifiability.expected_informative_reads` rather than a second copy of the
    formula.  An earlier version reimplemented it and hardcoded a mean effective length
    of 1000 against the package's 1500, which made the simulation draw 1.5x the fragments
    the prediction assumed -- a 0.82x error in the realised SD, comfortably inside the
    checkpoint's own tolerance and therefore silent.  One formula, one place.

    Flat in TPM within a class, as the survey's own precision figures assume; the two
    classes of ``primary_comparison`` are scaled to meet ``ratio`` while the gene's
    expected total is held fixed, and everything else in the gene keeps the base rate.
    """
    seqs, groups, pc = gene["seqs"], gene["groups"], gene["pc"]
    base = {t: I.expected_informative_reads(
        1.0, tpm, len(s), depth=30e6 * depth_scale,
        mean_efflen=mean_efflen, frag_mean=FRAG_MEAN) for t, s in seqs.items()}

    in_a = set(groups[pc[0]])
    in_b = set(groups[pc[1]])
    mass_a = sum(base[t] for t in in_a)
    mass_b = sum(base[t] for t in in_b)
    if ratio is not None and mass_a > 0 and mass_b > 0:
        # hold (mass_a + mass_b) fixed, move the split to `ratio`
        total = mass_a + mass_b
        want_a = total * ratio / (1.0 + ratio)
        want_b = total - want_a
        for t in in_a:
            base[t] *= want_a / mass_a
        for t in in_b:
            base[t] *= want_b / mass_b

    if expectations:
        return base
    return {t: int(rng.poisson(v)) for t, v in base.items()}


def draw_fragments(seq, n, rng, bias_5p=0.0):
    """Fragment lengths and start positions for n fragments from one transcript.

    ``bias_5p`` > 0 pulls starts toward the 5' end; < 0 is the exact mirror toward the 3'
    end, which is the direction degraded or polyA-selected RNA actually skews; 0 is
    uniform.  Off by default so its effect stays separable.

    The 3' case reflects the integer start rather than the uniform draw.  The obvious
    ``u = 1.0 - u`` looks equivalent and is not: ``1.0 - u`` rounds to exactly 1.0 for any
    ``u < 2**-53``, which at |bias| = 2 happens for one draw in 2e5 -- about 19
    out-of-range starts per run, measured, not estimated.  Reflecting the integer cannot
    leave [0, span-1] at all, and consumes the same rng either way so the fragment counts
    stay identical across bias values and the runs remain paired.
    """
    if n <= 0:
        return np.empty(0, int), np.empty(0, int)
    length = rng.normal(FRAG_MEAN, FRAG_SD, n).astype(int)
    np.clip(length, READ_LEN, max(READ_LEN, len(seq)), out=length)
    span = np.maximum(1, len(seq) - length + 1)
    if bias_5p:
        u = rng.random(n) ** (1.0 + abs(bias_5p))
        start = (u * span).astype(int)
        if bias_5p < 0:
            start = span - 1 - start
    else:
        start = (rng.random(n) * span).astype(int)
    return start, length


def mutate(reads, rng, rate=None):
    """Substitution error, vectorised over an (n, read_len) character array.

    ``rate`` is read from the module at call time rather than bound as a default, so that
    setting ``ERROR_RATE`` actually turns error off -- a default argument captures the
    value at import and silently ignores a later change.
    """
    rate = ERROR_RATE if rate is None else rate
    if rate <= 0:
        return reads
    hit = rng.random(reads.shape) < rate
    if hit.any():
        reads[hit] = BASES[rng.integers(0, 4, hit.sum())]
    return reads


def emit(seqs, counts, rng, f1, f2, tag, bias_5p=0.0):
    """Returns the per-transcript counts actually written, which is what the truth records.

    Fragments can be dropped -- a transcript shorter than a read, a fragment whose ends do
    not both fit -- and a truth file that recorded what was *drawn* would then credit the
    quantifier with reads that were never in the library.  None of the 1911 transcripts in
    the current survey is short enough for this to fire, which is exactly why it has to be
    handled here rather than remembered.
    """
    written = {}
    qual = QUAL * READ_LEN
    for tid, n in counts.items():
        written[tid] = 0
        if n <= 0:
            continue
        seq = seqs[tid]
        if len(seq) < READ_LEN:
            continue
        start, length = draw_fragments(seq, n, rng, bias_5p)
        frag_end = np.minimum(start + length, len(seq))
        arr = np.frombuffer(seq.encode(), dtype="S1")
        for i in range(n):
            a, b = int(start[i]), int(frag_end[i])
            if b - a < READ_LEN:
                a = max(0, b - READ_LEN)
            r1 = arr[a:a + READ_LEN].astype("U1").copy()
            r2raw = arr[max(0, b - READ_LEN):b].astype("U1").copy()
            if r1.size < READ_LEN or r2raw.size < READ_LEN:
                continue
            r1 = mutate(r1.reshape(1, -1), rng)[0]
            r2raw = mutate(r2raw.reshape(1, -1), rng)[0]
            s1 = "".join(r1)
            s2 = revcomp("".join(r2raw))
            name = "%s:%s:%d" % (tag, tid, i)
            f1.write("@%s/1\n%s\n+\n%s\n" % (name, s1, qual))
            f2.write("@%s/2\n%s\n+\n%s\n" % (name, s2, qual))
            written[tid] += 1
    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir")
    ap.add_argument("outdir")
    ap.add_argument("--depth-scale", type=float, default=1.0)
    ap.add_argument("--ratio", type=float, default=None,
                    help="true class A : class B, with the gene total held fixed. Omit "
                         "for the neutral case -- flat TPM, which is what the package's "
                         "own precision figure assumes. Setting it moves the Poisson "
                         "weights, so the prediction must be recomputed from these "
                         "abundances rather than taken from analyze()")
    ap.add_argument("--replicate", type=int, default=0)
    ap.add_argument("--tpm", type=float, default=10.0)
    ap.add_argument("--bias-5p", type=float, default=0.0,
                    help="positional coverage skew: >0 toward the 5' end, <0 toward "
                         "the 3' end (degraded / polyA-selected RNA), 0 uniform")
    ap.add_argument("--mean-efflen", type=float, default=I.DEFAULT_MEAN_EFFLEN,
                    help="must match what the prediction is computed with; the default "
                         "is the package's own")
    ap.add_argument("--genes", default=None)
    a = ap.parse_args()

    only = set(a.genes.split(",")) if a.genes else None
    genes = load_genes(a.rundir, only)
    if not genes:
        sys.exit("no genes found in %s" % a.rundir)
    os.makedirs(a.outdir, exist_ok=True)
    rng = np.random.default_rng(a.replicate * 7919 + 13)

    truth = {"replicate": a.replicate, "depth_scale": a.depth_scale, "ratio": a.ratio,
             "tpm": a.tpm, "bias_5p": a.bias_5p, "read_length": READ_LEN,
             "frag_mean": FRAG_MEAN, "frag_sd": FRAG_SD, "error_rate": ERROR_RATE,
             "mean_efflen": a.mean_efflen,
             "genes": {}}
    base = os.path.join(a.outdir, "rep%d" % a.replicate)
    total = 0
    with gzip.open(base + "_1.fq.gz", "wt") as f1, gzip.open(base + "_2.fq.gz", "wt") as f2:
        for gene, g in sorted(genes.items()):
            counts = transcript_fragment_counts(g, a.tpm, a.depth_scale, a.ratio, rng,
                                                a.mean_efflen)
            wrote = emit(g["seqs"], counts, rng, f1, f2,
                         "%s_r%d" % (gene, a.replicate), a.bias_5p)
            n = sum(wrote.values())
            total += n
            lost = sum(counts.values()) - n
            cls = {c: sum(wrote.get(t, 0) for t in g["groups"][c]) for c in g["pc"]}
            truth["genes"][gene] = {
                "primary_comparison": g["pc"],
                "fragments_per_transcript": wrote,
                "fragments_drawn_not_written": lost,
                "class_fragments": cls,
                "true_log2_ratio": (math.log2(cls[g["pc"][0]] / cls[g["pc"][1]])
                                    if cls[g["pc"][0]] > 0 and cls[g["pc"][1]] > 0
                                    else None),
                "fragments_written": n,
            }
    truth["fragments_written"] = total
    json.dump(truth, open(base + ".truth.json", "w"))
    print("replicate %d: %d genes, %d fragment pairs -> %s_{1,2}.fq.gz"
          % (a.replicate, len(genes), total, base))


if __name__ == "__main__":
    main()
