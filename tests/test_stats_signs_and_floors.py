"""Stouffer's sign and weights, p underflow, every floor, the stratified variance, and the
fold interval at n <= 5 (issue #20).

Each test names the number the audit computed on 2.5.0 and the number the corrected rule
gives; both were recomputed here before the code was changed.
"""
import csv
import itertools
import json

import numpy as np
import pytest
import scipy.stats
from scipy.stats import norm

from isoform_dominance import _selftest, cli, extract, stats

CFG = {"gene": "G", "groups": {"short": ["T1"], "long": ["T2"]},
       "primary_comparison": ["short", "long"]}
SIX_UP = (np.arange(6) + 10.0, np.ones(6))             # 6/6, exact p 0.03125


def _perdonor(path, A, B):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cohort", "donor", "condition", "short_TPM", "long_TPM"])
        for i, (a, b) in enumerate(zip(A, B, strict=True)):
            w.writerow(["C", "d%d" % i, "control", a, b])
    return str(path)


def _run(tmp_path, cohorts, **kw):
    paths = {name: _perdonor(tmp_path / ("%s.csv" % name), A, B)
             for name, (A, B) in cohorts.items()}
    kw.setdefault("n_boot", 0)
    return stats.run(CFG, "control", paths, str(tmp_path / "out"), **kw)


def _rows(tmp_path):
    with open(tmp_path / "out_stats.csv", newline="") as f:
        return list(csv.DictReader(f))


# --------------------------------------------------------------------------- #
# A1: the sign is the direction of the signed-rank statistic
# --------------------------------------------------------------------------- #
def test_a_median_below_zero_is_a_negative_sign():
    """d = [1, 1, -5, -6]: median -2, W+ 3 of 10.  The count of positive pairs is two of
    four, a tie, which 2.5.0 gave to +."""
    assert stats.signed_rank_direction([1.0, 1.0, -5.0, -6.0]) == -1.0
    det = stats.paired_stat_detail([2.0, 2.0, 1.0, 1.0], [1.0, 1.0, 6.0, 7.0], n_boot=0)
    assert det["n_greater"] * 2 >= det["n"]                # the old rule's +
    assert det["direction"] == -1.0


def test_counting_and_ranking_part_and_the_combination_follows_the_ranks(tmp_path):
    """A = [2,2,2,1,1,1], B = [1,1,1,21,21,21]: three pairs up by 1, three down by 20, median
    -9.5, W+ 6 of 21.  Combined with a 6/6 cohort, 2.5.0's count sign gave P 0.0283; the
    direction of the statistic that produced the cohort's p gives 0.393."""
    one = (np.array([2.0, 2, 2, 1, 1, 1]), np.array([1.0, 1, 1, 21, 21, 21]))
    res = _run(tmp_path, {"mixed": one, "up": SIX_UP})
    mixed = res["detail"][0]
    assert np.median(one[0] - one[1]) == -9.5
    assert mixed["direction"] == -1.0
    assert res["combination"]["stouffer"]["p"] == pytest.approx(0.39335791897434136, rel=1e-9)
    assert res["combination"]["stouffer"]["p"] != pytest.approx(0.0283, abs=1e-3)


def test_two_up_three_zero_none_down_is_positive():
    """2.5.0 put n_greater * 2 >= n over all five pairs, zeros included, and gave -."""
    det = stats.paired_stat_detail([2.0, 3.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0],
                                   n_boot=0)
    assert det["n_greater"] * 2 < det["n"]
    assert det["direction"] == 1.0


def test_a_balanced_statistic_has_no_direction():
    assert stats.signed_rank_direction([1.0, -1.0, 2.0, -2.0]) == 0.0
    assert stats.signed_rank_direction([0.0, 0.0]) == 0.0


def _scipy_direction(d, zero_method):
    """+1 when SciPy's one-sided 'greater' p is the smaller, -1 for 'less', 0 for neither."""
    ps = {alt: scipy.stats.wilcoxon(d, zero_method=zero_method, alternative=alt).pvalue
          for alt in ("greater", "less")}
    return float(np.sign(ps["less"] - ps["greater"]))


@pytest.mark.parametrize("zero_method", ["wilcox", "pratt", "zsplit"])
def test_the_direction_is_the_one_scipys_own_statistic_takes(zero_method):
    """Checked against SciPy's one-sided tests, not against a re-derivation: with zeros,
    pratt and zsplit rank them, which can move the direction away from wilcox's."""
    rng = np.random.default_rng(4)
    seen = set()
    for _ in range(80):
        n = int(rng.integers(2, 9))
        d = rng.choice([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0], size=n)
        if not np.any(d):
            continue
        got = stats.signed_rank_direction(d, zero_method)
        assert got == _scipy_direction(d, zero_method), (d, zero_method)
        seen.add(got)
    assert seen == {-1.0, 0.0, 1.0}


def test_pratt_can_point_where_wilcox_has_no_direction():
    d = [0.0, 0.0, 0.0, 1.0, 1.0, -3.0]
    assert stats.signed_rank_direction(d, "wilcox") == 0.0
    assert stats.signed_rank_direction(d, "pratt") == 1.0


# --------------------------------------------------------------------------- #
# A2: the weight is the number of pairs the test used
# --------------------------------------------------------------------------- #
def test_zero_differences_carry_no_weight(tmp_path):
    """Two informative pairs, both up, and eight zeros (p 0.5), with a 6/6 cohort: weighted by
    sqrt(n) 2.5.0 gave P 0.064, by sqrt(n_effective) 0.0276."""
    ten = (np.array([2.0, 3.0] + [1.0] * 8), np.ones(10))
    res = _run(tmp_path, {"zeros": ten, "up": SIX_UP})
    assert res["detail"][0]["n_effective"] == 2
    assert res["detail"][0]["p"] == 0.5
    assert res["combination"]["stouffer"]["p"] == pytest.approx(0.02762612316161699, rel=1e-9)


def test_pratt_keeps_the_zeros_in_the_weight(tmp_path):
    ten = (np.array([2.0, 3.0] + [1.0] * 8), np.ones(10))
    res = _run(tmp_path, {"zeros": ten, "up": SIX_UP}, zero_method="pratt")
    det = res["detail"][0]
    assert det["n_effective"] == 10
    z = (np.sqrt(10) * norm.isf(det["p"] / 2) + np.sqrt(6) * norm.isf(0.03125 / 2)) / 4.0
    assert res["combination"]["stouffer"]["z"] == pytest.approx(z, rel=1e-12)


# --------------------------------------------------------------------------- #
# A3: a p of zero
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("p", [0.0, 5e-324, 1e-320])
def test_a_p_of_zero_gives_a_finite_z(p):
    """SciPy's exact signed-rank p is 0.0 from about 55 untied pairs (54 on SciPy 1.17),
    and p / 2 underflows to 0 for the smallest subnormal; norm.isf(0) is inf."""
    up, down = stats._signed_z(p, 1.0), stats._signed_z(p, -1.0)
    assert np.isfinite(up) and up > 30
    assert down == -up
    one = stats.stouffer([("a", 6, 1.0, p)])
    assert np.isfinite(one["z"]) and one["p"] > 0.0


def test_two_opposed_cohorts_at_p_zero_cancel_instead_of_giving_nan():
    got = stats.stouffer([("a", 6, 1.0, 0.0), ("b", 6, -1.0, 0.0)])
    assert abs(got["z"]) < 1e-9                  # was inf - inf
    assert got["p"] == pytest.approx(1.0)


def test_a_missing_p_is_not_read_as_the_most_extreme_one():
    assert np.isnan(stats._signed_z(float("nan"), 1.0))


# --------------------------------------------------------------------------- #
# A4: the floor counts only the differences whose sign can flip
# --------------------------------------------------------------------------- #
def _ranks(d, zero_method):
    """The ranks each zero method gives |d|, and the differences they belong to."""
    if zero_method == "wilcox":
        d = d[d != 0]
    return d, scipy.stats.rankdata(np.abs(d))


def _enumerated_minimum(d, zero_method):
    """Smallest two-sided p the exact sign-flip null of ``zero_method`` can return for |d|.

    Every sign pattern of the non-zero differences is equally likely; a zero has no sign to
    flip under any zero method.  zsplit gives each side half a zero's rank, a constant."""
    d, r = _ranks(np.asarray(d, dtype=float), zero_method)
    nz = np.flatnonzero(d != 0)
    stats_ = []
    for signs in itertools.product((1.0, -1.0), repeat=nz.size):
        s = np.zeros(d.size)
        s[nz] = signs
        stats_.append(r[s > 0].sum() + (0.5 * r[s == 0].sum() if zero_method == "zsplit"
                                        else 0.0))
    w = np.sort(np.array(stats_))
    at_most = np.searchsorted(w, w, side="right") / w.size        # P(W <= x)
    at_least = 1.0 - np.searchsorted(w, w, side="left") / w.size  # P(W >= x)
    return float(np.min(np.minimum(1.0, 2.0 * np.minimum(at_most, at_least))))


@pytest.mark.parametrize("zero_method", ["wilcox", "pratt", "zsplit"])
def test_the_floor_is_the_enumerated_minimum_at_every_n_up_to_ten(zero_method):
    """n = 6 with one zero under pratt: 2.5.0 reported 0.03125, the minimum over the 32 sign
    patterns is 0.0625 = 2^(1 - n_nonzero)."""
    checked = 0
    for n in range(1, 11):
        for zeros in range(0, min(n, 4)):
            for tied in (False, True):
                mags = np.ones(n - zeros) * 2.0 if tied else np.arange(1.0, n - zeros + 1)
                d = np.concatenate([np.zeros(zeros), mags])
                det = stats.paired_stat_detail(d + 5.0, np.full(n, 5.0), n_boot=0,
                                               zero_method=zero_method)
                assert det["p_floor"] == pytest.approx(_enumerated_minimum(d, zero_method),
                                                       rel=1e-12), (n, zeros, tied)
                assert det["n_nonzero"] == n - zeros
                checked += 1
    assert checked == 68


def test_the_pratt_case_of_the_audit():
    A, B = np.array([1.0, 2, 3, 4, 5, 6]), np.ones(6)
    det = stats.paired_stat_detail(A, B, n_boot=0, zero_method="pratt")
    assert det["n_effective"] == 6 and det["n_nonzero"] == 5
    assert det["p_floor"] == 0.0625
    assert det["p"] >= det["p_floor"]


# --------------------------------------------------------------------------- #
# A5: the stratified statistic's variance under ties, and its floor
# --------------------------------------------------------------------------- #
def test_tied_absolute_differences_shrink_the_variance():
    """Five equal |d|: 13.75 without the correction, 11.25 with it (t = 5: 120/48 = 2.5)."""
    got = stats.stratified_signed_rank([np.full(5, 2.0)])
    assert got["z"] == pytest.approx(7.5 / np.sqrt(11.25), rel=1e-12)
    assert got["z"] != pytest.approx(7.5 / np.sqrt(13.75), rel=1e-6)


@pytest.mark.parametrize("d", [
    [2.0, 2.0, 2.0, 2.0, 2.0],
    [1.0, 1.0, 2.0, 3.0, 3.0, 3.0, -4.0, 5.0],
    [0.5, -0.5, 1.5, 1.5, 2.5, -2.5, 2.5, 3.0, 0.0],
])
def test_one_stratum_is_scipys_tie_corrected_normal_approximation(d):
    """With one stratum the weight cancels; SciPy's asymptotic p, without continuity
    correction, carries the tie-corrected variance."""
    ref = scipy.stats.wilcoxon(d, zero_method="wilcox", method="asymptotic",
                               correction=False)
    assert stats.stratified_signed_rank([np.array(d)])["p"] == pytest.approx(ref.pvalue,
                                                                             rel=1e-9)


def test_an_approximate_p_below_the_attainable_minimum_is_flagged():
    """d = [1..5]: P 0.0431 by the normal approximation, while every sign pattern of five
    pairs gives at least 0.0625."""
    got = stats.stratified_signed_rank([np.arange(1.0, 6.0)])
    assert got["p"] == pytest.approx(0.04311444678307529, rel=1e-9)
    assert got["p_floor"] == 0.0625
    assert got["underpowered"] is True
    assert got["approximation_below_floor"] is True
    assert got["note"] == stats.BELOW_FLOOR_NOTE == \
        "approximation below the attainable exact minimum"


def test_the_stratified_floor_counts_the_untied_pairs_of_every_stratum():
    got = stats.stratified_signed_rank([np.array([1.0, -2.0, 0.0]),
                                        np.array([3.0, 4.0, 5.0])])
    assert got["p_floor"] == 2.0 ** (1 - 5)
    assert got["p"] > got["p_floor"]
    assert got["approximation_below_floor"] is False
    assert got["note"] is None


def test_the_flag_reaches_the_csv_and_the_json(tmp_path, capsys):
    res = _run(tmp_path, {"a": (np.arange(1.0, 6.0) + 1, np.ones(5))})
    assert res["combination"]["stratified_signed_rank"]["approximation_below_floor"] is True
    row = {r["cohort"]: r for r in _rows(tmp_path)}["COMBINED_STRATIFIED_SIGNED_RANK"]
    assert row["note"] == stats.BELOW_FLOOR_NOTE
    assert row["resolution_floor_P"] == "0.0625"


def test_a_per_cohort_approximation_below_its_floor_is_flagged_too():
    """pratt at n = 14 with ten zeros: the normal approximation SciPy uses there reads 0.0465,
    below the 0.125 four flippable signs can reach."""
    A, B = np.concatenate([np.ones(10), [2.0, 3.0, 4.0, 5.0]]), np.ones(14)
    det = stats.paired_stat_detail(A, B, n_boot=0, zero_method="pratt")
    assert det["wilcoxon_method"] == "asymptotic"
    assert det["p"] < det["p_floor"] == 0.125
    assert det["approximation_below_floor"] is True
    clean = stats.paired_stat_detail(*SIX_UP, n_boot=0)
    assert clean["approximation_below_floor"] is False


# --------------------------------------------------------------------------- #
# A6: the combination rows carry their floor
# --------------------------------------------------------------------------- #
def _selftest_run(tmp_path, n_boot=0):
    info = _selftest.generate(str(tmp_path / "work"))
    perdonor = {}
    for cohort, p in info.items():
        out = str(tmp_path / ("pd_%s.csv" % cohort))
        extract.run(_selftest.CONFIG, p["quantdir"], p["samplemap"], cohort, out)
        perdonor[cohort] = out
    return stats.run(_selftest.CONFIG, "control", perdonor, str(tmp_path / "out"),
                     n_boot=n_boot)


def test_every_row_of_the_table_has_a_floor(tmp_path):
    """Both self-test cohorts sit at their floor (5/5, 6/6), so the Stouffer P is its own
    floor; the stratified floor is 2^(1-11)."""
    res = _selftest_run(tmp_path)
    rows = {r["cohort"]: r for r in _rows(tmp_path)}
    assert all(r["resolution_floor_P"] for r in rows.values())
    stouffer = res["combination"]["stouffer"]
    assert stouffer["p_floor"] == stouffer["p"]
    assert rows["COMBINED_STOUFFER"]["resolution_floor_P"] == "0.004419" \
        == rows["COMBINED_STOUFFER"]["paired_wilcoxon_P"]
    assert rows["COMBINED_STRATIFIED_SIGNED_RANK"]["resolution_floor_P"] == "0.0009766"
    assert rows["COMBINED_STOUFFER"]["underpowered"] == "no"
    assert rows["COMBINED_STRATIFIED_SIGNED_RANK"]["underpowered"] == "no"


def test_the_stouffer_floor_is_every_cohort_at_its_floor_in_one_direction():
    per = [("a", 5, -1.0, 0.5, 0.0625), ("b", 6, 1.0, 0.03125, 0.03125)]
    z_max = (np.sqrt(5) * norm.isf(0.0625 / 2) + np.sqrt(6) * norm.isf(0.03125 / 2)) \
        / np.sqrt(11)
    got = stats.stouffer(per)
    assert got["p_floor"] == pytest.approx(2 * norm.sf(z_max), rel=1e-12)
    assert got["p"] > got["p_floor"]
    assert got["underpowered"] is False


def test_three_donors_a_cohort_cannot_reach_005_however_combined(tmp_path):
    """Two cohorts of three: each floor 0.25, Z_max = 1.627, floor 0.104."""
    three = (np.array([4.0, 5.0, 6.0]), np.ones(3))
    res = _run(tmp_path, {"a": three, "b": three})
    stouffer = res["combination"]["stouffer"]
    assert stouffer["p_floor"] == pytest.approx(0.1037, abs=1e-4)
    assert stouffer["underpowered"] is True
    assert {r["cohort"]: r for r in _rows(tmp_path)}["COMBINED_STOUFFER"]["underpowered"] \
        == "yes"


def test_a_stouffer_without_every_floor_reports_none():
    assert np.isnan(stats.stouffer([("a", 6, 1.0, 0.03125)])["p_floor"])


# --------------------------------------------------------------------------- #
# A7: at n <= 5 the bootstrap interval of the median is the sample range
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_the_bootstrap_interval_of_a_median_of_five_or_fewer_is_the_range(n):
    """Why the label: P(bootstrap median = min) is 0.25, 7/27, 13/256 and 181/3125 at n = 2..5,
    every one above 0.025, so the 95% percentile interval is the sample's min and max."""
    x = np.array([4.93, 22.34, 25.94, 28.26, 32.76])[:n]
    lo, hi = stats._bootstrap_median_ci(x, 10_000, 0, 0.95)
    assert (lo, hi) == (x.min(), x.max())


def test_at_six_the_bootstrap_interval_is_inside_the_range():
    x = np.array([29.6, 31.0, 36.0, 39.4, 47.0, 54.32])
    lo, hi = stats._bootstrap_median_ci(x, 10_000, 0, 0.95)
    assert x.min() < lo and hi < x.max()


def test_gse228458_reports_its_range_and_says_so(tmp_path):
    """The bundled (4.93, 32.76) is the min and max of the five donor ratios."""
    res = _selftest_run(tmp_path, n_boot=2000)
    by = {d["cohort"]: d for d in res["detail"]}
    small, six = by["GSE228458"], by["GSE137619"]
    with open(tmp_path / "pd_GSE228458.csv") as f:
        ratios = [float(r["short_TPM"]) / float(r["long_TPM"]) for r in csv.DictReader(f)]
    assert small["fold_ci"] == (min(ratios), max(ratios))
    assert small["fold_ci_method"] == "range (n<=5)"
    assert six["fold_ci_method"] == "bootstrap"
    assert res["pooled"]["fold_ci_method"] == "bootstrap"
    rows = {r["cohort"]: r for r in _rows(tmp_path)}
    assert (rows["GSE228458"]["fold_CI_low"], rows["GSE228458"]["fold_CI_high"]) == \
        ("4.93", "32.76")
    assert rows["GSE228458"]["fold_CI_method"] == "range (n<=5)"
    assert rows["GSE137619"]["fold_CI_method"] == "bootstrap"
    assert rows["COMBINED_STOUFFER"]["fold_CI_method"] == ""


def test_no_interval_without_replicates_and_no_method_either():
    det = stats.paired_stat_detail([2.0, 4.0, 3.0], [1.0, 2.0, 1.0], n_boot=0)
    assert np.isnan(det["fold_ci"][0]) and det["fold_ci_method"] is None


def test_the_new_columns_are_appended_and_the_old_ones_keep_their_place(tmp_path):
    _selftest_run(tmp_path)
    with open(tmp_path / "out_stats.csv", newline="") as f:
        head = next(csv.reader(f))
    assert head == ["cohort", "n", "short>long", "median_fold", "fold_CI_low",
                    "fold_CI_high", "paired_wilcoxon_P", "resolution_floor_P",
                    "underpowered", "wilcoxon_method", "fold_CI_method", "note"]


# --------------------------------------------------------------------------- #
# what the reader sees
# --------------------------------------------------------------------------- #
def test_the_stats_command_prints_each_combinations_floor(tmp_path, capsys):
    paths = {"a": _perdonor(tmp_path / "a.csv", np.arange(1.0, 6.0) + 1, np.ones(5))}
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(CFG))
    assert cli.main(["stats", "--config", str(cfg), "--perdonor", "a=%s" % paths["a"],
                     "--out", str(tmp_path / "out")]) == 0
    out = capsys.readouterr().out
    line = {ln.split()[0]: ln for ln in out.splitlines() if ln.strip()}
    assert "floor 0.0625" in line["STOUFFER"]
    assert "floor 0.0625" in line["STRATIFIED_SIGNED_RANK"]
    assert stats.BELOW_FLOOR_NOTE in line["STRATIFIED_SIGNED_RANK"]
    assert "range" in line["a"]


def test_the_selftest_shows_the_floor_of_each_combination(capsys):
    assert cli.main(["selftest"]) == 0
    lines = capsys.readouterr().out.splitlines()
    stouffer = next(ln for ln in lines if ln.strip().startswith("Stouffer"))
    assert "P=0.004419" in stouffer and "floor=0.004419" in stouffer
    assert cli.main(["selftest", "--json"]) == 0
    got = json.loads(capsys.readouterr().out)["combinations"]
    assert got["stouffer"]["p_floor"] == got["stouffer"]["p"]
    assert got["stratified_signed_rank"]["p_floor"] == 2.0 ** -10
