"""Paired isoform-group statistics + figure.

The estimand is a within-donor contrast: each donor contributes one long-class and
one short-class abundance, and the question is whether the class ranking is
consistent across donors.  Three things about that are easy to get wrong and are
handled explicitly here.

*Combining cohorts.*  Pooling donors from independent cohorts into one signed-rank
test treats cohort as if it did not exist, when in practice cohorts differ in tissue
handling, library preparation and depth.  The pooled test is still reported, because
it is what earlier releases reported and what the bundled reference result is stated
in, but two stratified combinations are reported beside it -- a weighted Stouffer
combination of the per-cohort *exact* tests, and a weighted combination of the
within-stratum signed-rank statistics -- and those are the numbers to quote when the
cohorts are genuinely independent studies.  Stouffer is the default of the two because
at these sample sizes the per-stratum exact p-value is trustworthy and the normal
approximation behind a combined rank statistic is not.

*The exact-test floor.*  A two-sided exact signed-rank test on ``n`` untied pairs
cannot return a p-value below ``2^(1-n)``: at n = 5 the floor is 0.0625, so no
arrangement of five donors is significant at 0.05.  Reporting a non-significant
p-value without that context invites the reader to conclude the effect is absent when
the design could never have shown it.  The floor is computed and reported alongside
every test, and flagged when it exceeds 0.05.

*Effect size.*  A median fold-change with no interval is a point estimate presented
as if it were a measurement; a donor-level bootstrap interval is reported with it.

Matplotlib is imported lazily inside ``run`` so that importing this module (and
therefore the CLI) does not pay the matplotlib/font-cache startup cost for
subcommands that never plot (``--version``, ``annotate``, ``identifiability``).
"""
import csv
import math
import os

import numpy as np
from scipy.stats import norm, wilcoxon

from .io import primary_pair

PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]

#: Bootstrap replicates for the fold-change interval, and the seed that makes the
#: interval reproducible across runs and machines.
DEFAULT_N_BOOT = 10_000
DEFAULT_SEED = 0

#: How cohorts are combined by default in the reported summary.
DEFAULT_COMBINATION = "stouffer"


def load_perdonor(path, condition, gA, gB):
    don, A, B = [], [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            if condition not in (None, "", "all") and r["condition"] != condition:
                continue
            don.append(r["donor"]); A.append(float(r["%s_TPM" % gA])); B.append(float(r["%s_TPM" % gB]))
    return don, np.array(A), np.array(B)


# --------------------------------------------------------------------------- #
# the exact-test floor
# --------------------------------------------------------------------------- #
def signed_rank_resolution_floor(n, alternative="two-sided"):
    """Finest p-value the signed-rank permutation null can resolve with ``n`` pairs.

    Under the exhaustive sign-permutation null each of the ``n`` non-zero differences
    takes either sign with probability 1/2, giving ``2^n`` equally likely assignments.
    Exactly one of them puts every difference on the same side, so a two-sided test
    cannot report anything below ``2 / 2^n = 2^(1-n)``.  At n = 5 that is 0.0625: five
    donors in perfect agreement still cannot reject at 0.05, and reporting the
    non-significant p-value without that context invites the reader to conclude the
    effect is absent when the design could never have shown it.

    What ``n`` must be
        The count of pairs with a **non-zero** difference.  ``zero_method="wilcox"``
        discards tied pairs, so the value to pass is the effective one --
        :func:`paired_stat_detail` computes it and reports it as ``n_effective``.

    What ties do
        Nothing, to this bound.  Ties among the absolute differences change the shape of
        the null distribution and therefore the attainable p-values *between* the
        extremes, but the extreme itself is still reached by exactly one sign
        assignment, because midranks preserve the rank sum.  (Checked against SciPy: a
        five-pair comparison whose differences are all equal and all positive still
        returns exactly 0.0625.)  An earlier version of this docstring claimed the floor
        rises under ties; it does not.

    What this is not
        It is a property of the *design*, not of the procedure that was run.  It says
        what an exhaustive sign permutation could resolve; it does not certify that the
        p-value beside it came from an exact computation.  ``paired_stat_detail``
        reports ``exact_null_clean`` for that.
    """
    if n <= 0:
        return float("nan")
    one_sided = 0.5 ** n
    return 2.0 * one_sided if alternative == "two-sided" else one_sided


# --------------------------------------------------------------------------- #
# per-cohort statistics
# --------------------------------------------------------------------------- #
def paired_stat(A, B):
    """Return (n, n_A>B, two-sided exact Wilcoxon P, median fold A/B).

    Edge cases are handled explicitly: an empty input returns NaNs; if every pair
    is tied (no nonzero differences) the signed-rank test is undefined and P is
    NaN; the median fold-change ignores non-finite ratios (e.g. zero denominators).

    Kept as a four-tuple for callers written against v2.1; :func:`paired_stat_detail`
    returns the same numbers plus the floor, the bootstrap interval and the tie count.
    """
    d = paired_stat_detail(A, B, n_boot=0)
    return d["n"], d["n_greater"], d["p"], d["median_fold"]


def paired_stat_detail(A, B, n_boot=DEFAULT_N_BOOT, seed=DEFAULT_SEED,
                       zero_method="wilcox", ci=0.95):
    """Paired comparison with its floor, effect-size interval and tie accounting.

    ``zero_method`` is passed through to :func:`scipy.stats.wilcoxon` and stated
    explicitly rather than left to the SciPy default, because the choice changes the
    effective ``n`` -- ``"wilcox"`` discards tied pairs -- and therefore changes the
    achievable floor.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    n = len(A)
    out = {"n": n, "n_greater": 0, "p": float("nan"), "median_fold": float("nan"),
           "fold_ci": (float("nan"), float("nan")), "n_ties": 0,
           "n_effective": 0, "p_floor": float("nan"),
           "underpowered": False, "zero_method": zero_method}
    if n == 0:
        return out

    d = A - B
    out["n_greater"] = int(np.sum(A > B))
    out["n_ties"] = int(np.sum(d == 0))
    n_eff = n - out["n_ties"] if zero_method == "wilcox" else n
    out["n_effective"] = int(n_eff)
    nz = d[d != 0]
    # ties among |d| leave the floor alone but do change the null's shape, and some
    # SciPy versions decline to compute an exact p in their presence; record it so a
    # reader is not left assuming the reported p is the permutation one
    out["exact_null_clean"] = bool(nz.size and
                                   np.unique(np.abs(nz)).size == nz.size)
    out["p_floor"] = signed_rank_resolution_floor(n_eff)
    out["underpowered"] = bool(n_eff > 0 and out["p_floor"] > 0.05)

    if np.allclose(A, B):
        out["p"] = float("nan")            # all differences zero -> test undefined
    else:
        try:
            out["p"] = float(wilcoxon(A, B, alternative="two-sided",
                                      zero_method=zero_method).pvalue)
        except ValueError:
            out["p"] = float("nan")

    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = A / B
    finite = ratios[np.isfinite(ratios)]
    out["median_fold"] = float(np.median(finite)) if finite.size else float("nan")

    if n_boot and finite.size > 1:
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, finite.size, size=(int(n_boot), finite.size))
        boots = np.median(finite[idx], axis=1)
        lo, hi = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
        out["fold_ci"] = (float(np.quantile(boots, lo)), float(np.quantile(boots, hi)))
    return out


# --------------------------------------------------------------------------- #
# combining cohorts
# --------------------------------------------------------------------------- #
def _signed_z(p, direction):
    """Two-sided p plus a direction -> a signed z-score."""
    if not np.isfinite(p) or p <= 0.0:
        p = np.nextafter(0, 1)
    p = min(p, 1.0)
    return float(direction) * float(norm.isf(p / 2.0))


def stouffer(per_cohort):
    """Weighted Stouffer combination of independent per-cohort exact tests.

    Each cohort contributes its own exact two-sided p-value, signed by the direction
    of its median difference and weighted by ``sqrt(n)``.  Combining p-values rather
    than pooling donors keeps each cohort's test exact -- which matters at these
    sample sizes, where the normal approximation used by a stratified rank statistic
    is poor -- while still refusing to treat donors from different studies as
    exchangeable.

    ``per_cohort`` is an iterable of ``(name, n, direction, p)``.
    """
    zs, ws = [], []
    for _, n, direction, p in per_cohort:
        if not n or not np.isfinite(p):
            continue
        zs.append(_signed_z(p, direction))
        ws.append(math.sqrt(n))
    if not zs:
        return {"z": float("nan"), "p": float("nan"), "k": 0}
    z = float(np.dot(ws, zs) / math.sqrt(float(np.dot(ws, ws))))
    return {"z": z, "p": float(2.0 * norm.sf(abs(z))), "k": len(zs)}


def stratified_signed_rank(strata):
    """Weighted stratified combination of within-stratum signed-rank statistics.

    Within each stratum the positive-rank sum ``W+`` is computed on the untied pairs and
    centred and scaled by its null moments; strata are then combined with weights
    ``1/(n+1)``.

    **On the name.**  Those weights are van Elteren's design-free choice, but van
    Elteren's test itself is a stratified *two-sample* (Wilcoxon rank-sum) procedure for
    independent groups within strata.  The data here are paired within donor, so the
    within-stratum statistic is the signed-rank one and this is *not* van Elteren's
    test; only the weighting scheme is borrowed, and the function is named for what it
    does rather than for him.

    A normal approximation is used for the combined statistic, so with very few donors
    per stratum this should be read alongside, not instead of, :func:`stouffer`, which
    keeps each stratum's test exact.

    ``strata`` is an iterable of paired difference arrays.
    """
    num = den = 0.0
    used = 0
    for d in strata:
        d = np.asarray(d, dtype=float)
        d = d[d != 0]
        n = d.size
        if n == 0:
            continue
        order = np.argsort(np.abs(d), kind="mergesort")
        ranks = np.empty(n, dtype=float)
        ranks[order] = np.arange(1, n + 1, dtype=float)
        # average ranks over ties in |d|
        absd = np.abs(d)
        for v in np.unique(absd):
            m = absd == v
            if m.sum() > 1:
                ranks[m] = ranks[m].mean()
        w_plus = float(ranks[d > 0].sum())
        mean = n * (n + 1) / 4.0
        var = n * (n + 1) * (2 * n + 1) / 24.0
        weight = 1.0 / (n + 1.0)
        num += weight * (w_plus - mean)
        den += weight ** 2 * var
        used += 1
    if used == 0 or den <= 0:
        return {"z": float("nan"), "p": float("nan"), "k": used}
    z = num / math.sqrt(den)
    return {"z": float(z), "p": float(2.0 * norm.sf(abs(z))), "k": used}


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def run(config, condition, cohorts, out, n_boot=DEFAULT_N_BOOT, seed=DEFAULT_SEED,
        zero_method="wilcox"):
    """cohorts: {name: perdonor.csv}. Writes <out>.{png,pdf,svg} + <out>_stats.csv.

    The returned dict keeps ``per_cohort`` and ``combined`` in their v2.1 shapes --
    ``combined`` is still the donor-pooled test, so a caller pinned to that number
    (the bundled self-test, and the reference result quoted in the documentation)
    sees no change -- and adds ``detail``, ``combination`` and ``p_floor``.
    """
    import matplotlib as mpl; mpl.use("Agg")
    import matplotlib.pyplot as plt
    gene = config.get("gene", "gene")
    gA, gB = primary_pair(config)
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 9, "pdf.fonttype": 42, "svg.fonttype": "none"})
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    names = list(cohorts)
    fig, axes = plt.subplots(1, len(names), figsize=(3.6 * len(names), 3.6), squeeze=False)
    allA, allB, statrows, details, diffs = [], [], [], [], []
    for i, name in enumerate(names):
        don, A, B = load_perdonor(cohorts[name], condition, gA, gB)
        if len(A) == 0:
            raise ValueError(
                "cohort %s: no donors matched condition %r in %s."
                % (name, condition, cohorts[name]))
        allA += list(A); allB += list(B)
        det = paired_stat_detail(A, B, n_boot=n_boot, seed=seed, zero_method=zero_method)
        det["cohort"] = name
        details.append(det)
        diffs.append(A - B)
        n, ngt, p, fold = det["n"], det["n_greater"], det["p"], det["median_fold"]
        statrows.append((name, n, ngt, p, fold))
        ax = axes[0, i]; c = PALETTE[i % len(PALETTE)]; fl = 1e-3
        for k in range(n):
            ax.plot([0, 1], [max(A[k], fl), max(B[k], fl)], color="#999", lw=0.8, zorder=1)
        ax.scatter(np.zeros(n), np.clip(A, fl, None), s=34, c=c, edgecolor="white", lw=0.5, zorder=3)
        ax.scatter(np.ones(n), np.clip(B, fl, None), s=34, c="#9aa0a6", edgecolor="white", lw=0.5, zorder=3)
        ax.set_yscale("log"); ax.set_xlim(-0.4, 1.4); ax.set_xticks([0, 1]); ax.set_xticklabels([gA, gB])
        ax.set_ylabel("%s TPM (log)" % gene); ax.set_title("%s (%s n=%d)" % (name, condition, n), fontsize=9.5)
        label = "%d/%d  %.0fx  P=%.3f" % (ngt, n, fold, p)
        if det["underpowered"]:
            label += "  (floor %.3f)" % det["p_floor"]
        ax.text(0.5, 1.16, label, transform=ax.transAxes,
                ha="center", va="top", fontsize=8, color="#333")
        ax.spines[["top", "right"]].set_visible(False)

    pooled = paired_stat_detail(np.array(allA), np.array(allB), n_boot=n_boot,
                               seed=seed, zero_method=zero_method)
    cn, cgt, cp, cfold = pooled["n"], pooled["n_greater"], pooled["p"], pooled["median_fold"]

    combo = {
        "stouffer": stouffer([(d["cohort"], d["n"],
                               1.0 if d["n_greater"] * 2 >= d["n"] else -1.0, d["p"])
                              for d in details]),
        "stratified_signed_rank": stratified_signed_rank(diffs),
        "pooled": {"z": float("nan"), "p": cp, "k": len(details)},
    }
    headline = combo.get(DEFAULT_COMBINATION, combo["pooled"])

    fig.suptitle("%s: %s vs %s - combined %s n=%d, %d/%d, P=%.4g (%s P=%.4g)"
                 % (gene, gA, gB, condition, cn, cgt, cn, cp,
                    DEFAULT_COMBINATION, headline["p"]),
                 fontsize=10.5, fontweight="bold", y=1.04)
    for ext in ("png", "pdf", "svg"):
        fig.savefig("%s.%s" % (out, ext), dpi=300, bbox_inches="tight")
    plt.close(fig)

    with open(out + "_stats.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cohort", "n", "%s>%s" % (gA, gB), "median_fold",
                    "fold_CI_low", "fold_CI_high", "paired_wilcoxon_P",
                    "resolution_floor_P", "underpowered"])
        for d in details:
            w.writerow([d["cohort"], d["n"], "%d/%d" % (d["n_greater"], d["n"]),
                        "%.2f" % d["median_fold"], "%.2f" % d["fold_ci"][0],
                        "%.2f" % d["fold_ci"][1], "%.4g" % d["p"],
                        "%.4g" % d["p_floor"], "yes" if d["underpowered"] else "no"])
        w.writerow(["COMBINED_POOLED", cn, "%d/%d" % (cgt, cn), "%.2f" % cfold,
                    "%.2f" % pooled["fold_ci"][0], "%.2f" % pooled["fold_ci"][1],
                    "%.4g" % cp, "%.4g" % pooled["p_floor"],
                    "yes" if pooled["underpowered"] else "no"])
        for key in ("stouffer", "stratified_signed_rank"):
            w.writerow(["COMBINED_%s" % key.upper(), cn, "", "", "", "",
                        "%.4g" % combo[key]["p"], "", ""])

    return {"per_cohort": statrows, "combined": (cn, cgt, cp, cfold),
            "detail": details, "pooled": pooled, "combination": combo,
            "headline_combination": DEFAULT_COMBINATION}
