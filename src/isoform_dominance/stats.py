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

*The exact-test floor.*  A two-sided exact signed-rank test on ``n`` pairs with a
non-zero difference cannot return a p-value below ``2^(1-n)``: at n = 5 the floor is
0.0625, so no arrangement of five donors is significant at 0.05.  Reporting a
non-significant p-value without that context invites the reader to conclude the effect
is absent when the design could never have shown it.  Every test is reported with its
floor -- each cohort's, the donor-pooled one, and the two combinations -- and flagged
when the floor exceeds 0.05.  A p-value from a normal approximation can fall below the
floor its own data could reach exactly; that is flagged as well (:data:`BELOW_FLOOR_NOTE`).

*Effect size.*  A median fold-change with no interval is a point estimate presented
as if it were a measurement; a donor-level bootstrap interval is reported with it.  With
five ratios or fewer the bootstrap percentile interval of the median is the sample's
minimum and maximum, so it is computed and labelled as that range.

Matplotlib is imported lazily inside ``run`` so that importing this module (and
therefore the CLI) does not pay the matplotlib/font-cache startup cost for
subcommands that never plot (``--version``, ``annotate``, ``identifiability``).
"""
import csv
import math
import os
import warnings

import numpy as np
import scipy.stats
from scipy.stats import norm, wilcoxon

from .io import InputError, primary_pair

#: Class colours. Identity follows the isoform class, not the panel index -- the
#: cohort is already encoded by which panel a donor is in, so reusing the colour
#: channel for cohort would leave the same class drawn in different hues from one
#: panel to the next. Validated as a categorical pair against a light surface
#: (CVD dE 24.7, normal-vision dE 33.6), and both classes are direct-labelled on the
#: x axis, so identity never rests on colour alone.
CLASS_COLORS = ("#2a78d6", "#eb6834")
INK = "#3d3d38"
INK_MUTED = "#84837b"
HAIRLINE = "#d4d3cc"

#: Retained for callers written against v2.1.
PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]

#: Bootstrap replicates for the fold-change interval, and the seed that makes the
#: interval reproducible across runs and machines.
DEFAULT_N_BOOT = 10_000
DEFAULT_SEED = 0

#: How cohorts are combined by default in the reported summary.  Stouffer's inputs are
#: the per-cohort exact tests, where the stratified signed-rank alternative reaches its P
#: through a normal approximation that is not validated at the stratum sizes this package
#: is built for; the alternative is reported alongside, never in place of, the default.
#: Neither is claimed to hold its level at every alpha: the inverse-normal transform of a
#: discrete exact p-value is not a normal variate.  Each reports its floor, the smallest P
#: it can return given the cohorts' sizes, so a combination that sits at its floor shows it.
DEFAULT_COMBINATION = "stouffer"

#: Below this many finite fold ratios the reported interval is the sample's range: the
#: 95% bootstrap percentile interval of a median of five or fewer values is its minimum
#: and maximum, because the bootstrap median equals the minimum with probability 0.25,
#: 7/27, 13/256 and 181/3125 at n = 2 to 5, each above 0.025 (0.0087 at n = 6).
RANGE_MAX_N = 5
FOLD_CI_BOOTSTRAP = "bootstrap"
FOLD_CI_RANGE = "range (n<=5)"

#: What a row says when its P comes from a normal approximation and lies below the floor:
#: no sign pattern of the same pairs could give that P exactly.
BELOW_FLOOR_NOTE = "approximation below the attainable exact minimum"

#: Font stack for the figures, first match wins.  Arial leads because it is what most
#: journals ask for; DejaVu Sans is Matplotlib's bundled fallback and is what a machine
#: without Arial will use.  Set this to a single family before calling :func:`run` to
#: make a figure reproduce identically across machines -- ``scripts/make_docs_example.py``
#: pins DejaVu Sans for exactly that reason.
FONT_STACK = ["Arial", "DejaVu Sans"]

#: Display names for the combination methods, used in the figure and the CLI.
COMBINATION_LABELS = {
    "stouffer": "Stouffer",
    "stratified_signed_rank": "stratified signed-rank",
    "pooled": "donor-pooled",
}


def load_perdonor(path, condition, gA, gB):
    don, A, B = [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        need = ["donor", "%s_TPM" % gA, "%s_TPM" % gB]
        if condition not in (None, "", "all"):
            need.insert(1, "condition")
        absent = [c for c in need if c not in (reader.fieldnames or [])]
        if absent:
            raise InputError("per-donor table %s has no %s column (its columns: %s); "
                             "`extract` writes one for each group of the config"
                             % (path, ", ".join(absent), ", ".join(reader.fieldnames or [])
                                or "none"))
        for r in reader:
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
        The count of pairs with a **non-zero** difference, under every zero method.
        ``zero_method="wilcox"`` discards the zero differences; ``"pratt"`` and
        ``"zsplit"`` rank them, but a zero has no sign to flip, so it adds no pattern to
        the null and the extreme is still one pattern in ``2^n``.  Through 2.6.0 the
        floor under ``pratt`` counted every pair: at n = 6 with one zero it read 0.03125
        where the 32 sign patterns reach no lower than 0.0625.
        :func:`paired_stat_detail` reports the count as ``n_nonzero``.

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
    """Return (n, n_A>B, two-sided Wilcoxon P, median fold A/B).

    P is SciPy's ``method="auto"`` one: exact only without zeros or ties, see
    :func:`paired_stat_detail`.

    Edge cases are handled explicitly: an empty input returns NaNs; if every pair
    is tied (no nonzero differences) the signed-rank test is undefined and P is
    NaN; the median fold-change ignores non-finite ratios (e.g. zero denominators).

    Kept as a four-tuple for callers written against v2.1; :func:`paired_stat_detail`
    returns the same numbers plus the floor, the bootstrap interval and the tie count.
    """
    d = paired_stat_detail(A, B, n_boot=0)
    return d["n"], d["n_greater"], d["p"], d["median_fold"]


def _wilcoxon_auto_method(d):
    """The computation ``wilcoxon(method="auto")`` selects for differences ``d``.

    SciPy does not report which one it used, so its dispatch (SciPy >= 1.15, the
    declared floor) is mirrored here: ``"asymptotic"`` when n > 50; otherwise
    ``"exact"`` when there are neither ties nor zeros, an exhaustive ``"permutation"``
    test when n <= 13, and ``"asymptotic"`` above that.  ``n`` counts every pair, zeros
    included, as SciPy's does.  ``"asymptotic"`` is the normal approximation, without
    continuity correction.

    Earlier SciPy used the normal approximation for any zero difference at any n, which
    is why the floor is 1.15.
    """
    d = np.asarray(d, dtype=float)
    n = d.size
    has_zeros = bool(np.any(d == 0))
    if n > 50:
        return "asymptotic"
    nz = np.abs(d[d != 0])
    has_ties = np.unique(nz).size != nz.size
    if not (has_ties or has_zeros):
        return "exact"
    return "permutation" if n <= 13 else "asymptotic"


def _wilcoxon_method_arg(label):
    """The explicit SciPy ``method=`` argument a method label stands for."""
    if label == "permutation":
        return scipy.stats.PermutationMethod()
    return label


def signed_rank_direction(d, zero_method="wilcox"):
    """Direction of the signed-rank statistic for differences ``d``: +1.0, -1.0 or 0.0.

    The sign of ``W+ - W-``, the rank sums of the positive and the negative differences,
    with the ranks the test itself uses: average ranks of ``|d|`` over the non-zero
    differences for ``zero_method="wilcox"``, over every difference for ``"pratt"`` and
    ``"zsplit"`` (a zero's rank goes to neither side, or half to each, which moves both
    sums alike).  Under the sign-flip null ``W+`` is centred on ``(W+ + W-)/2``, so this is
    the side of its centre the statistic that produced the p-value fell on -- for
    ``wilcox``, the sign of ``W+ - n(n+1)/4`` over the ``n`` non-zero pairs.  0.0 when it
    sits on the centre, where the test states no direction.

    Not the sign of the median difference, nor of the count of positive pairs, which
    2.6.0 used: ``d = [1, 1, -5, -6]`` has two of four pairs positive and a median of -2,
    and its ``W+`` is 3 of 10.
    """
    d = np.asarray(d, dtype=float)
    if zero_method == "wilcox":
        d = d[d != 0]
    if not np.any(d):
        return 0.0
    ranks = scipy.stats.rankdata(np.abs(d))
    diff = float(ranks[d > 0].sum() - ranks[d < 0].sum())
    return float(np.sign(diff))


def _bootstrap_median_ci(finite, n_boot, seed, ci):
    """The donor-bootstrap percentile interval of the median of ``finite``."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, finite.size, size=(int(n_boot), finite.size))
    boots = np.median(finite[idx], axis=1)
    lo, hi = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    return float(np.quantile(boots, lo)), float(np.quantile(boots, hi))


def _below_floor(p, floor):
    """Is ``p`` below ``floor``, beyond rounding?  False when either is not a number."""
    return bool(np.isfinite(p) and np.isfinite(floor) and p < floor * (1.0 - 1e-9))


def paired_stat_detail(A, B, n_boot=DEFAULT_N_BOOT, seed=DEFAULT_SEED,
                       zero_method="wilcox", ci=0.95):
    """Paired comparison with its floor, effect-size interval and tie accounting.

    ``zero_method`` is passed through to :func:`scipy.stats.wilcoxon` and stated
    explicitly rather than left to the SciPy default, because the choice changes the
    effective ``n`` -- ``"wilcox"`` discards tied pairs -- and therefore changes the
    achievable floor.

    The p-value is SciPy's ``method="auto"`` one, unchanged, and ``auto`` does not
    always mean exact: with a zero difference or a tie it runs an exhaustive
    permutation test at n <= 13 and the normal approximation above that.
    ``wilcoxon_method`` names the computation that
    actually produced ``p`` -- ``"exact"``, ``"permutation"`` or ``"asymptotic"``.  The
    name comes from :func:`_wilcoxon_auto_method` and is reported only if re-running
    SciPy with that method named explicitly reproduces ``p`` bit for bit; otherwise it
    is ``"unresolved"``, so a SciPy release that changes its rule shows up as that
    rather than as a wrong label.  It is ``None`` when no test was run.

    ``n_effective`` is the number of pairs the test ranks -- the non-zero ones under
    ``wilcox``, all of them under ``pratt`` and ``zsplit`` -- and weights the cohort in
    :func:`stouffer`.  ``n_nonzero`` is the number whose sign can flip, and sets
    ``p_floor`` under every zero method.  ``direction`` is
    :func:`signed_rank_direction`'s.  ``approximation_below_floor`` is true when ``p``
    lies below ``p_floor``, which only a normal approximation can do.

    ``fold_ci`` is the donor-bootstrap 95% interval of the median ratio, with
    ``fold_ci_method`` ``"bootstrap"``; from five finite ratios or fewer it is their
    minimum and maximum, ``"range (n<=5)"``, which is what the bootstrap percentile
    interval of such a median is at ``ci >= 0.95`` (:data:`RANGE_MAX_N`).  Neither is
    computed with ``n_boot=0`` or fewer than two finite ratios, and ``fold_ci_method`` is
    then None.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    n = len(A)
    out = {"n": n, "n_greater": 0, "p": float("nan"), "median_fold": float("nan"),
           "fold_ci": (float("nan"), float("nan")), "fold_ci_method": None, "n_ties": 0,
           "n_effective": 0, "n_nonzero": 0, "direction": 0.0, "p_floor": float("nan"),
           "underpowered": False, "approximation_below_floor": False,
           "zero_method": zero_method, "wilcoxon_method": None}
    if n == 0:
        return out

    d = A - B
    out["n_greater"] = int(np.sum(A > B))
    out["n_ties"] = int(np.sum(d == 0))
    n_eff = n - out["n_ties"] if zero_method == "wilcox" else n
    out["n_effective"] = int(n_eff)
    out["n_nonzero"] = int(n - out["n_ties"])
    out["direction"] = signed_rank_direction(d, zero_method)
    nz = d[d != 0]
    # ties among |d| leave the floor alone but do change the null's shape; record it so
    # a reader is not left assuming the reported p is the permutation one.
    # `wilcoxon_method` below says which computation SciPy actually ran.
    out["exact_null_clean"] = bool(nz.size and
                                   np.unique(np.abs(nz)).size == nz.size)
    out["p_floor"] = signed_rank_resolution_floor(out["n_nonzero"])
    out["underpowered"] = bool(out["n_nonzero"] > 0 and out["p_floor"] > 0.05)

    if np.allclose(A, B):
        out["p"] = float("nan")            # all differences zero -> test undefined
    else:
        try:
            out["p"] = float(wilcoxon(A, B, alternative="two-sided",
                                      zero_method=zero_method, method="auto").pvalue)
        except ValueError:
            out["p"] = float("nan")
        if np.isfinite(out["p"]):
            label = _wilcoxon_auto_method(d)
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")   # already raised by the call above
                    again = float(wilcoxon(A, B, alternative="two-sided",
                                           zero_method=zero_method,
                                           method=_wilcoxon_method_arg(label)).pvalue)
            except ValueError:
                again = float("nan")
            out["wilcoxon_method"] = label if again == out["p"] else "unresolved"
    out["approximation_below_floor"] = _below_floor(out["p"], out["p_floor"])

    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = A / B
    finite = ratios[np.isfinite(ratios)]
    out["median_fold"] = float(np.median(finite)) if finite.size else float("nan")

    if n_boot and finite.size > 1:
        if finite.size <= RANGE_MAX_N and ci >= 0.95:
            # the bootstrap interval, computed without the bootstrap: see RANGE_MAX_N
            out["fold_ci"] = (float(finite.min()), float(finite.max()))
            out["fold_ci_method"] = FOLD_CI_RANGE
        else:
            out["fold_ci"] = _bootstrap_median_ci(finite, n_boot, seed, ci)
            out["fold_ci_method"] = FOLD_CI_BOOTSTRAP
    return out


# --------------------------------------------------------------------------- #
# combining cohorts
# --------------------------------------------------------------------------- #
def _signed_z(p, direction):
    """Two-sided p plus a direction -> a signed z-score, always finite for a p in [0, 1].

    ``p / 2`` is held at the smallest normal double (2.2e-308, z = 37.5): a p of 0.0 --
    SciPy's exact signed-rank p is 0.0 from about 55 untied pairs -- or one whose half
    underflows gave ``norm.isf(0) = inf`` through 2.6.0, and an infinite z swamps every
    other cohort, or turns two opposed cohorts into NaN.  A p that is not a number gives
    NaN, not the most extreme z."""
    if not np.isfinite(p):
        return float("nan")
    p = min(max(float(p), 0.0), 1.0)
    return float(direction) * float(norm.isf(max(p / 2.0, _TINY)))


#: The smallest normal double, below which ``p / 2`` is not taken (see :func:`_signed_z`).
_TINY = float(np.finfo(float).tiny)


def stouffer(per_cohort):
    """Weighted Stouffer combination of independent per-cohort exact tests.

    Each cohort contributes its own two-sided p-value as a z-score, signed by the
    direction of the signed-rank statistic that produced it
    (:func:`signed_rank_direction`, ``direction`` in :func:`paired_stat_detail`) and
    weighted by ``sqrt(n)``, ``n`` being the number of pairs that test ranked
    (``n_effective``: a pair with no difference carries no weight under ``wilcox``).
    Combining p-values rather than pooling donors keeps each cohort's test exact --
    which matters at these sample sizes, where the normal approximation used by a
    stratified rank statistic is poor -- while still refusing to treat donors from
    different studies as exchangeable.

    Through 2.6.0 the sign was that of the count of positive pairs, ties going to +, and
    the weight ``sqrt`` of every donor: with one cohort of three pairs up by 1 and three
    down by 20 beside a 6/6 cohort, P read 0.0283 where the cohorts' own statistics give
    0.393; a cohort of two pairs up and eight with no difference weighed as ten.

    ``per_cohort`` is an iterable of ``(name, n, direction, p)`` or ``(name, n,
    direction, p, p_floor)``.  With every cohort's floor given, ``p_floor`` is the
    combination's: each cohort at its floor in one direction, ``Z_max = sum(w_i z_i) /
    sqrt(sum(w_i^2))`` with ``z_i`` the z of the cohort's floor, and ``P = 2 sf(Z_max)``;
    ``underpowered`` when that exceeds 0.05.  Otherwise ``p_floor`` is NaN.
    """
    zs, ws, tops = [], [], []
    for _, n, direction, p, *floor in per_cohort:
        if not n or not np.isfinite(p):
            continue
        zs.append(_signed_z(p, direction))
        ws.append(math.sqrt(n))
        tops.append(_signed_z(floor[0], 1.0) if floor else float("nan"))
    if not zs:
        return {"z": float("nan"), "p": float("nan"), "k": 0, "p_floor": float("nan"),
                "underpowered": False, "approximation_below_floor": False}
    norm_w = math.sqrt(float(np.dot(ws, ws)))
    z = float(np.dot(ws, zs) / norm_w)
    p = float(2.0 * norm.sf(abs(z)))
    floor = float(2.0 * norm.sf(np.dot(ws, tops) / norm_w))
    return {"z": z, "p": p, "k": len(zs), "p_floor": floor,
            "underpowered": bool(np.isfinite(floor) and floor > 0.05),
            "approximation_below_floor": _below_floor(p, floor)}


def stratified_signed_rank(strata):
    """Weighted stratified combination of within-stratum signed-rank statistics.

    Within each stratum the positive-rank sum ``W+`` is computed on the untied pairs and
    centred and scaled by its null moments; strata are then combined with weights
    ``1/(n+1)``.  The null variance of ``W+`` with average ranks over tied ``|d|`` is
    ``n(n+1)(2n+1)/24 - sum(t^3 - t)/48`` over the tie groups, ``t`` the size of each;
    through 2.6.0 the second term was left out, which overstated the variance and so
    understated ``|z|`` whenever a stratum had ties.

    **On the name.**  Those weights are van Elteren's design-free choice, but van
    Elteren's test itself is a stratified *two-sample* (Wilcoxon rank-sum) procedure for
    independent groups within strata.  The data here are paired within donor, so the
    within-stratum statistic is the signed-rank one and this is *not* van Elteren's
    test; only the weighting scheme is borrowed, and the function is named for what it
    does rather than for him.

    A normal approximation is used for the combined statistic, so with very few donors
    per stratum this should be read alongside, not instead of, :func:`stouffer`, which
    keeps each stratum's test exact.  ``p_floor`` is what the exact stratified sign-flip
    null could return at its extreme -- every pair of every stratum on one side, one
    pattern in ``2^N`` -- so ``2^(1 - N)`` with ``N`` the untied pairs of all strata;
    ``underpowered`` when it exceeds 0.05.  The approximation can fall below it: one
    stratum of five pairs all up gives 0.0431 against a floor of 0.0625.
    ``approximation_below_floor`` says so, and ``note`` is then :data:`BELOW_FLOOR_NOTE`.

    ``strata`` is an iterable of paired difference arrays.
    """
    num = den = 0.0
    used = n_total = 0
    for d in strata:
        d = np.asarray(d, dtype=float)
        d = d[d != 0]
        n = d.size
        if n == 0:
            continue
        absd = np.abs(d)
        ranks = scipy.stats.rankdata(absd)          # average ranks over ties in |d|
        w_plus = float(ranks[d > 0].sum())
        mean = n * (n + 1) / 4.0
        _, t = np.unique(absd, return_counts=True)
        var = n * (n + 1) * (2 * n + 1) / 24.0 - float(np.sum(t ** 3 - t)) / 48.0
        weight = 1.0 / (n + 1.0)
        num += weight * (w_plus - mean)
        den += weight ** 2 * var
        used += 1
        n_total += n
    floor = signed_rank_resolution_floor(n_total)
    out = {"z": float("nan"), "p": float("nan"), "k": used, "p_floor": floor,
           "underpowered": bool(n_total > 0 and floor > 0.05),
           "approximation_below_floor": False, "note": None}
    if used == 0 or den <= 0:
        return out
    z = num / math.sqrt(den)
    out.update(z=float(z), p=float(2.0 * norm.sf(abs(z))))
    if _below_floor(out["p"], floor):
        out.update(approximation_below_floor=True, note=BELOW_FLOOR_NOTE)
    return out


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def run(config, condition, cohorts, out, n_boot=DEFAULT_N_BOOT, seed=DEFAULT_SEED,
        zero_method="wilcox"):
    """cohorts: {name: perdonor.csv}. Writes <out>.{png,pdf,svg} + <out>_stats.csv.

    The returned dict keeps ``per_cohort`` and ``combined`` in their v2.1 shapes --
    ``combined`` is still the donor-pooled test, so a caller pinned to that number
    (the bundled self-test, and the reference result quoted in the documentation)
    sees no change -- and adds ``detail`` (each cohort's :func:`paired_stat_detail`),
    ``pooled`` (the same for every donor), ``combination`` and
    ``headline_combination``.  Each entry of ``combination`` has ``z``, ``p``, ``k``,
    ``p_floor``, ``underpowered`` and ``approximation_below_floor``; the stratified one
    has ``note`` too.  ``z`` of ``pooled`` is NaN: that test is no z-combination.

    Every row of ``<out>_stats.csv`` carries its floor in ``resolution_floor_P``, the
    two combination rows included, and ``underpowered`` when it exceeds 0.05.  Two
    columns follow the v2.6 ones: ``fold_CI_method`` (``bootstrap``, ``range (n<=5)``,
    or empty with no interval) and ``note`` (:data:`BELOW_FLOOR_NOTE` when a P from a
    normal approximation lies below its floor).
    """
    import matplotlib as mpl; mpl.use("Agg")
    import matplotlib.pyplot as plt
    gene = config.get("gene", "gene")
    gA, gB = primary_pair(config)
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": list(FONT_STACK),
                         "font.size": 9, "pdf.fonttype": 42, "svg.fonttype": "none"})
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    names = list(cohorts)
    # shared y across panels: these are small multiples of the same measurement, and
    # an independent scale per panel makes the slopes look comparable when they are not
    fig, axes = plt.subplots(1, len(names), figsize=(3.6 * len(names), 4.0),
                             squeeze=False, sharey=True)
    fig.subplots_adjust(top=0.80, bottom=0.16, wspace=0.10)
    allA, allB, statrows, details, diffs = [], [], [], [], []
    for i, name in enumerate(names):
        don, A, B = load_perdonor(cohorts[name], condition, gA, gB)
        if len(A) == 0:
            raise InputError(
                "cohort %s: no donors matched condition %r in %s."
                % (name, condition, cohorts[name]))
        allA += list(A); allB += list(B)
        det = paired_stat_detail(A, B, n_boot=n_boot, seed=seed, zero_method=zero_method)
        det["cohort"] = name
        details.append(det)
        diffs.append(A - B)
        n, ngt, p, fold = det["n"], det["n_greater"], det["p"], det["median_fold"]
        statrows.append((name, n, ngt, p, fold))

        ax = axes[0, i]
        fl = 1e-3
        for k in range(n):
            ax.plot([0, 1], [max(A[k], fl), max(B[k], fl)],
                    color=HAIRLINE, lw=1.0, zorder=1, solid_capstyle="round")
        ax.scatter(np.zeros(n), np.clip(A, fl, None), s=46, c=CLASS_COLORS[0],
                   edgecolor="white", lw=1.2, zorder=3)
        ax.scatter(np.ones(n), np.clip(B, fl, None), s=46, c=CLASS_COLORS[1],
                   edgecolor="white", lw=1.2, zorder=3)

        ax.set_yscale("log")
        ax.set_xlim(-0.45, 1.45)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([gA, gB], fontsize=9.5, color=INK)
        ax.tick_params(axis="both", which="major", colors=INK_MUTED, length=3, width=0.8)
        ax.tick_params(axis="both", which="minor", colors=INK_MUTED, length=1.8, width=0.6)
        if i > 0:
            # sharey already suppresses the labels; the orphaned tick marks are noise.
            # which="both" matters here: a log axis carries minor ticks as well.
            ax.tick_params(axis="y", which="both", left=False, right=False)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(HAIRLINE)
            ax.spines[side].set_linewidth(0.8)
        ax.grid(axis="y", color=HAIRLINE, lw=0.6, alpha=0.6)
        ax.set_axisbelow(True)
        if i == 0:
            ax.set_ylabel("%s TPM  (log scale)" % gene, fontsize=9, color=INK)

        # the per-cohort numbers belong to the panel, so they live in its title --
        # floating them above the axes collided with the figure title
        sub = "%d/%d donors  ·  %.0fx" % (ngt, n, fold)
        lo, hi = det["fold_ci"]
        if np.isfinite(lo) and np.isfinite(hi):
            sub += " [%.0f-%.0f]" % (lo, hi)
        sub += "  ·  P = %.3g" % p
        ax.set_title("%s   (%s, n = %d)" % (name, condition, n),
                     fontsize=10, color=INK, pad=20, loc="center")
        ax.annotate(sub, xy=(0.5, 1.012), xycoords="axes fraction",
                    ha="center", va="bottom", fontsize=8.3, color=INK_MUTED)
        if det["underpowered"]:
            ax.annotate("the exact test cannot go below P = %.3g at n = %d"
                        % (det["p_floor"], det["n_nonzero"]),
                        xy=(0.5, -0.115), xycoords="axes fraction", ha="center",
                        va="top", fontsize=7.6, color=INK_MUTED, style="italic")

    pooled = paired_stat_detail(np.array(allA), np.array(allB), n_boot=n_boot,
                               seed=seed, zero_method=zero_method)
    cn, cgt, cp, cfold = pooled["n"], pooled["n_greater"], pooled["p"], pooled["median_fold"]

    combo = {
        "stouffer": stouffer([(d["cohort"], d["n_effective"], d["direction"], d["p"],
                               d["p_floor"]) for d in details]),
        "stratified_signed_rank": stratified_signed_rank(diffs),
        "pooled": {"z": float("nan"), "p": cp, "k": len(details),
                   "p_floor": pooled["p_floor"], "underpowered": pooled["underpowered"],
                   "approximation_below_floor": pooled["approximation_below_floor"]},
    }
    headline_key = DEFAULT_COMBINATION if DEFAULT_COMBINATION in combo else "pooled"
    fig.suptitle("%s: %s vs %s  ·  %d cohorts, %s donors"
                 % (gene, gA, gB, len(names), condition),
                 fontsize=11.5, fontweight="bold", color=INK, y=0.985)
    # one headline number in the header; the alternatives are named in a footnote so
    # they are available without three P values competing for the same slot
    fig.text(0.5, 0.925,
             "combined n = %d, %d/%d concordant  ·  %s P = %.3g"
             % (cn, cgt, cn, COMBINATION_LABELS.get(headline_key, headline_key),
                combo[headline_key]["p"]),
             ha="center", va="center", fontsize=8.6, color=INK_MUTED)
    alts = ["%s P = %.3g" % (COMBINATION_LABELS.get(k, k), combo[k]["p"])
            for k in ("stratified_signed_rank", "stouffer")
            if k != headline_key and np.isfinite(combo[k]["p"])]
    alts.append("donor-pooled P = %.3g" % cp)
    fig.text(0.995, 0.005, "also: " + "  ·  ".join(alts), ha="right", va="bottom",
             fontsize=7.4, color=INK_MUTED)

    for ext in ("png", "pdf", "svg"):
        fig.savefig("%s.%s" % (out, ext), dpi=300, bbox_inches="tight")
    plt.close(fig)

    with open(out + "_stats.csv", "w", newline="") as f:
        w = csv.writer(f)
        # new columns are appended last so that a reader indexing the older columns by
        # position is unaffected; the combination rows are not Wilcoxon tests
        w.writerow(["cohort", "n", "%s>%s" % (gA, gB), "median_fold",
                    "fold_CI_low", "fold_CI_high", "paired_wilcoxon_P",
                    "resolution_floor_P", "underpowered", "wilcoxon_method",
                    "fold_CI_method", "note"])
        for name, d, n_, ngt_ in ([(d["cohort"], d, d["n"], d["n_greater"]) for d in details]
                                  + [("COMBINED_POOLED", pooled, cn, cgt)]):
            w.writerow([name, n_, "%d/%d" % (ngt_, n_),
                        "%.2f" % d["median_fold"], "%.2f" % d["fold_ci"][0],
                        "%.2f" % d["fold_ci"][1], "%.4g" % d["p"],
                        "%.4g" % d["p_floor"], "yes" if d["underpowered"] else "no",
                        d["wilcoxon_method"] or "", d["fold_ci_method"] or "",
                        BELOW_FLOOR_NOTE if d["approximation_below_floor"] else ""])
        for key in ("stouffer", "stratified_signed_rank"):
            c = combo[key]
            w.writerow(["COMBINED_%s" % key.upper(), cn, "", "", "", "",
                        "%.4g" % c["p"], "%.4g" % c["p_floor"],
                        "yes" if c["underpowered"] else "no", "", "",
                        BELOW_FLOOR_NOTE if c["approximation_below_floor"] else ""])

    return {"per_cohort": statrows, "combined": (cn, cgt, cp, cfold),
            "detail": details, "pooled": pooled, "combination": combo,
            "headline_combination": DEFAULT_COMBINATION}
