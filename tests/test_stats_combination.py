"""The v2.2 statistics layer: exact-test floor, effect-size interval, cohort combination."""
import csv
import json

import numpy as np
import pytest
import scipy
import scipy.stats

from isoform_dominance import cli, stats


# --------------------------------------------------------------------------- #
# the exact-test floor
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("n,expected", [(5, 0.0625), (6, 0.03125), (11, 2 ** -10)])
def test_resolution_floor_matches_the_documented_value(n, expected):
    assert stats.signed_rank_resolution_floor(n) == pytest.approx(expected)


def test_resolution_floor_one_sided_is_half_the_two_sided():
    assert stats.signed_rank_resolution_floor(7, "greater") == pytest.approx(
        stats.signed_rank_resolution_floor(7) / 2)


def test_ties_do_not_raise_the_resolution_floor():
    """The floor is 2^(1-n) whether or not the absolute differences are tied.

    Only one sign assignment puts every difference on the same side, and midranks
    preserve the rank sum, so the extreme is reached exactly once either way. An earlier
    version of the docstring claimed ties push the floor up; this pins the correction.
    """
    tied = np.array([3.0, 4.0, 5.0, 6.0, 7.0]), np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    det = stats.paired_stat_detail(*tied)          # every difference is exactly 2.0
    assert det["n_effective"] == 5
    assert det["exact_null_clean"] is False        # |d| fully tied
    assert det["p"] == pytest.approx(2 ** -4)      # still the floor, not above it
    assert det["p_floor"] == pytest.approx(2 ** -4)


def test_zeros_shrink_the_effective_n_and_raise_the_floor():
    A = np.array([1.0, 4.0, 5.0, 6.0, 7.0])
    B = np.array([1.0, 2.0, 3.0, 4.0, 5.0])        # one tied pair -> dropped
    det = stats.paired_stat_detail(A, B)
    assert det["n_effective"] == 4
    assert det["p_floor"] == pytest.approx(2 ** -3)
    assert det["p"] == pytest.approx(2 ** -3)


def test_five_donors_are_flagged_underpowered_even_when_unanimous():
    A = np.array([10.0, 8.0, 6.0, 9.0, 7.0])
    B = np.array([1.0, 2.0, 1.5, 3.0, 0.5])
    det = stats.paired_stat_detail(A, B)
    assert det["n_greater"] == 5                     # every donor agrees
    assert det["p"] > 0.05                           # and it still cannot reach 0.05
    assert det["underpowered"] is True
    assert det["p_floor"] == pytest.approx(0.0625)


def test_eleven_donors_are_not_flagged_underpowered():
    rng = np.random.default_rng(0)
    A = rng.uniform(5, 10, 11)
    B = rng.uniform(0.1, 1.0, 11)
    det = stats.paired_stat_detail(A, B)
    assert det["underpowered"] is False
    assert det["p"] == pytest.approx(2 ** -10)


# --------------------------------------------------------------------------- #
# effect size
# --------------------------------------------------------------------------- #
def test_bootstrap_interval_brackets_the_point_estimate_and_is_deterministic():
    rng = np.random.default_rng(1)
    A = rng.uniform(4, 6, 12)
    B = rng.uniform(1, 2, 12)
    first = stats.paired_stat_detail(A, B, n_boot=2000, seed=7)
    second = stats.paired_stat_detail(A, B, n_boot=2000, seed=7)
    lo, hi = first["fold_ci"]
    assert lo <= first["median_fold"] <= hi
    assert first["fold_ci"] == second["fold_ci"]     # seeded, so reproducible


def test_bootstrap_can_be_disabled():
    det = stats.paired_stat_detail([2.0, 4.0], [1.0, 2.0], n_boot=0)
    assert np.isnan(det["fold_ci"][0])


def test_ties_are_counted_and_shrink_the_effective_n():
    A = np.array([1.0, 2.0, 3.0, 4.0])
    B = np.array([1.0, 1.0, 1.0, 1.0])
    det = stats.paired_stat_detail(A, B)
    assert det["n_ties"] == 1
    assert det["n_effective"] == 3                   # zero_method="wilcox" drops it
    assert det["p_floor"] == pytest.approx(0.25)


def test_paired_stat_tuple_api_is_unchanged():
    A, B = [3.0, 4.0, 5.0], [1.0, 1.0, 1.0]
    n, ngt, p, fold = stats.paired_stat(A, B)
    det = stats.paired_stat_detail(A, B)
    assert (n, ngt, fold) == (det["n"], det["n_greater"], det["median_fold"])
    assert p == pytest.approx(det["p"], nan_ok=True)


# --------------------------------------------------------------------------- #
# combining cohorts
# --------------------------------------------------------------------------- #
def test_stouffer_combines_concordant_cohorts_below_either_alone():
    per = [("c1", 5, 1.0, 0.0625), ("c2", 6, 1.0, 0.03125)]
    got = stats.stouffer(per)
    assert got["k"] == 2
    assert got["p"] < min(p for *_, p in per)


def test_stouffer_cancels_opposing_directions():
    per = [("c1", 6, 1.0, 0.03), ("c2", 6, -1.0, 0.03)]
    got = stats.stouffer(per)
    assert abs(got["z"]) < 1e-9
    assert got["p"] == pytest.approx(1.0)


def test_stouffer_ignores_cohorts_without_a_usable_p():
    per = [("c1", 6, 1.0, 0.03), ("empty", 0, 1.0, float("nan"))]
    assert stats.stouffer(per)["k"] == 1
    assert np.isnan(stats.stouffer([("x", 0, 1.0, float("nan"))])["p"])


def test_stratified_signed_rank_detects_a_consistent_within_stratum_shift():
    strata = [np.array([2.0, 3.0, 1.0, 4.0, 2.5]),
              np.array([1.0, 2.0, 3.0, 1.5, 2.0, 2.2])]
    got = stats.stratified_signed_rank(strata)
    assert got["k"] == 2
    assert got["z"] > 0
    assert got["p"] < 0.05


def test_stratified_signed_rank_is_null_for_symmetric_differences():
    strata = [np.array([1.0, -1.0, 2.0, -2.0])]
    got = stats.stratified_signed_rank(strata)
    assert got["p"] > 0.5


def test_stratified_signed_rank_skips_all_tied_strata():
    assert stats.stratified_signed_rank([np.zeros(5)])["k"] == 0
    assert np.isnan(stats.stratified_signed_rank([np.zeros(5)])["p"])


def test_stratified_signed_rank_is_invariant_to_rescaling_one_stratum():
    """Ranks are taken within a stratum, so a cohort's units cannot leak into another's.

    This is the property that makes a stratified statistic the right default here:
    pooling donors ranks a low-depth cohort's differences against a high-depth
    cohort's on a common scale, so a cohort with systematically larger TPMs takes the
    top ranks and dominates the result for reasons that have nothing to do with
    isoform usage.
    """
    rng = np.random.default_rng(3)
    s1 = np.array([0.03, 0.02, 0.05, 0.01, 0.04])       # unanimous, tiny magnitudes
    s2 = rng.normal(0, 1, 12)

    stratified = stats.stratified_signed_rank([s1, s2])
    rescaled = stats.stratified_signed_rank([s1 * 1000.0, s2])
    assert stratified["z"] == pytest.approx(rescaled["z"])

    pooled = stats.stratified_signed_rank([np.concatenate([s1, s2])])
    pooled_rescaled = stats.stratified_signed_rank([np.concatenate([s1 * 1000.0, s2])])
    assert pooled["z"] != pytest.approx(pooled_rescaled["z"])


def test_stratified_signed_rank_keeps_a_small_unanimous_stratum_visible():
    """A small cohort with tiny but perfectly consistent differences still counts."""
    rng = np.random.default_rng(3)
    small = np.array([0.03, 0.02, 0.05, 0.01, 0.04])     # unanimous, tiny
    large = rng.normal(0, 5, 30)                          # null, large magnitudes
    stratified = stats.stratified_signed_rank([small, large])
    pooled = stats.stratified_signed_rank([np.concatenate([small, large])])
    assert stratified["z"] > pooled["z"]


# --------------------------------------------------------------------------- #
# the driver
# --------------------------------------------------------------------------- #
def _perdonor(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cohort", "donor", "condition", "short_TPM", "long_TPM"])
        for donor, a, b in rows:
            w.writerow(["C", donor, "control", a, b])


def test_run_reports_pooled_and_stratified_side_by_side(tmp_path):
    p1, p2 = tmp_path / "c1.csv", tmp_path / "c2.csv"
    _perdonor(p1, [("d%d" % i, 5.0 + i, 1.0) for i in range(5)])
    _perdonor(p2, [("e%d" % i, 4.0 + i, 1.2) for i in range(6)])
    cfg = {"gene": "G", "groups": {"short": ["T1"], "long": ["T2"]},
           "primary_comparison": ["short", "long"]}
    res = stats.run(cfg, "control", {"c1": str(p1), "c2": str(p2)},
                    str(tmp_path / "out"), n_boot=500)

    # v2.1 shape preserved: `combined` is still the donor-pooled test
    assert res["combined"][0] == 11
    assert res["combined"][2] == pytest.approx(2 ** -10)
    # and the stratified alternatives are reported beside it
    for key in ("stouffer", "stratified_signed_rank", "pooled"):
        assert key in res["combination"]
    assert res["combination"]["stouffer"]["k"] == 2
    assert res["headline_combination"] == stats.DEFAULT_COMBINATION

    body = (tmp_path / "out_stats.csv").read_text()
    assert "resolution_floor_P" in body
    assert "COMBINED_POOLED" in body
    assert "COMBINED_STOUFFER" in body
    assert "COMBINED_STRATIFIED_SIGNED_RANK" in body
    assert "fold_CI_low" in body


# --------------------------------------------------------------------------- #
# which Wilcoxon method produced the p-value
# --------------------------------------------------------------------------- #
# SciPy's method="auto" picks the p-value computation from the data, and the rule
# changed in 1.15: before, ties were ignored and only zeros forced the normal
# approximation; from 1.15, ties or zeros select an exhaustive permutation test at
# n <= 13 and the normal approximation above that.
SCIPY_1_15 = tuple(int(x) for x in scipy.__version__.split(".")[:2]) >= (1, 15)

# 14 pairs, one zero difference, the rest distinct: every supported SciPy falls back to
# the normal approximation here, and the approximation is visibly not the exact value
ZERO14_B = np.full(14, 10.0)
ZERO14_A = ZERO14_B + np.array([0.0, 1.1, 2.3, -0.4, 3.7, 4.2, -1.9,
                                5.5, 6.1, 2.8, 7.3, -0.9, 8.6, 3.1])
# 6 pairs, no zeros, distinct |d|: exact on every supported SciPy
CLEAN6_B = np.full(6, 5.0)
CLEAN6_A = CLEAN6_B + np.array([1.5, 2.5, -0.5, 3.5, 4.5, 5.5])
# 6 pairs, no zeros, tied |d| with mixed signs: exact before 1.15, permutation after
TIED6_A = np.array([3.0, 4.0, 1.0, 6.0, 7.0, 8.0])
TIED6_B = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])


def _explicit(label):
    """The SciPy `method=` argument that a reported label stands for."""
    if label == "permutation":
        return scipy.stats.PermutationMethod()
    if label == "asymptotic" and not SCIPY_1_15:
        return "approx"
    return label


@pytest.mark.parametrize("A,B,expected", [
    (ZERO14_A, ZERO14_B, "asymptotic"),
    (CLEAN6_A, CLEAN6_B, "exact"),
    (TIED6_A, TIED6_B, "permutation" if SCIPY_1_15 else "exact"),
])
def test_wilcoxon_method_names_the_computation_that_produced_p(A, B, expected):
    det = stats.paired_stat_detail(A, B, n_boot=0)
    assert det["wilcoxon_method"] == expected
    # the value is still SciPy's own "auto" p, untouched ...
    auto = scipy.stats.wilcoxon(A, B, zero_method="wilcox", method="auto").pvalue
    assert det["p"] == auto
    # ... and the label is the explicit method that reproduces it
    named = scipy.stats.wilcoxon(A, B, zero_method="wilcox",
                                 method=_explicit(expected)).pvalue
    assert det["p"] == named


def _enumerated_exact_p(A, B):
    """Two-sided signed-rank p from every sign assignment of the non-zero differences.

    Computed here rather than asked of SciPy, because SciPy < 1.15 answers
    method="exact" with the normal approximation whenever a difference is zero.
    """
    d = np.asarray(A, dtype=float) - np.asarray(B, dtype=float)
    d = d[d != 0]
    ranks = scipy.stats.rankdata(np.abs(d))
    assert np.unique(ranks).size == ranks.size          # untied, so ranks are 1..n
    obs = ranks[d > 0].sum()
    signs = (np.arange(2 ** d.size)[:, None] >> np.arange(d.size)) & 1
    w = signs @ ranks
    return min(1.0, 2.0 * min(np.mean(w <= obs), np.mean(w >= obs)))


def test_enumerated_exact_p_agrees_with_scipy_where_scipy_is_exact():
    exact = scipy.stats.wilcoxon(CLEAN6_A, CLEAN6_B, method="exact").pvalue
    assert _enumerated_exact_p(CLEAN6_A, CLEAN6_B) == pytest.approx(exact, rel=1e-12)


def test_the_zero_difference_case_is_not_the_exact_p():
    """Why the label matters: here the reported p is not the exact one, silently."""
    det = stats.paired_stat_detail(ZERO14_A, ZERO14_B, n_boot=0)
    exact = _enumerated_exact_p(ZERO14_A, ZERO14_B)
    assert det["n_ties"] == 1
    assert det["p"] != pytest.approx(exact, rel=1e-3)
    assert det["wilcoxon_method"] != "exact"


def test_an_unconfirmed_method_label_is_not_reported(monkeypatch):
    """If SciPy's rule drifts from the one mirrored here, say so rather than guess."""
    monkeypatch.setattr(stats, "_wilcoxon_auto_method", lambda *a, **k: "exact")
    det = stats.paired_stat_detail(ZERO14_A, ZERO14_B, n_boot=0)
    assert det["wilcoxon_method"] == "unresolved"
    assert det["p"] == scipy.stats.wilcoxon(ZERO14_A, ZERO14_B).pvalue


def test_undefined_test_has_no_method():
    det = stats.paired_stat_detail([1.0, 2.0], [1.0, 2.0], n_boot=0)
    assert np.isnan(det["p"])
    assert det["wilcoxon_method"] is None


def _two_cohort_run(tmp_path):
    p1, p2 = tmp_path / "zero14.csv", tmp_path / "clean6.csv"
    _perdonor(p1, [("z%d" % i, a, b)
                   for i, (a, b) in enumerate(zip(ZERO14_A, ZERO14_B, strict=True))])
    _perdonor(p2, [("c%d" % i, a, b)
                   for i, (a, b) in enumerate(zip(CLEAN6_A, CLEAN6_B, strict=True))])
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"gene": "G", "groups": {"short": ["T1"], "long": ["T2"]},
                               "primary_comparison": ["short", "long"]}))
    return str(cfg), {"zero14": str(p1), "clean6": str(p2)}


def test_stats_csv_has_a_wilcoxon_method_column(tmp_path):
    cfg, cohorts = _two_cohort_run(tmp_path)
    stats.run(json.loads(open(cfg).read()), "control", cohorts,
              str(tmp_path / "out"), n_boot=0)
    with open(tmp_path / "out_stats.csv") as f:
        rows = {r["cohort"]: r for r in csv.DictReader(f)}
    assert rows["zero14"]["wilcoxon_method"] == "asymptotic"
    assert rows["clean6"]["wilcoxon_method"] == "exact"
    # 20 pooled pairs with one zero: the approximation on every supported SciPy
    assert rows["COMBINED_POOLED"]["wilcoxon_method"] == "asymptotic"


def test_stats_cli_prints_the_wilcoxon_method(tmp_path, capsys):
    cfg, cohorts = _two_cohort_run(tmp_path)
    rc = cli.main(["stats", "--config", cfg, "--n-boot", "0",
                   "--perdonor", "zero14=%s" % cohorts["zero14"],
                   "--perdonor", "clean6=%s" % cohorts["clean6"],
                   "--out", str(tmp_path / "out")])
    assert rc == 0
    lines = capsys.readouterr().out.splitlines()
    line = {ln.split()[0]: ln for ln in lines if ln.strip()}
    assert "(asymptotic)" in line["zero14"]
    assert "(exact)" in line["clean6"]
    assert "(asymptotic)" in line["POOLED"]
