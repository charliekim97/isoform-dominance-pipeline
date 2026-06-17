"""Edge cases for paired_stat, config validation, and extract input validation."""
import math
import numpy as np
import pytest
from isoform_dominance import stats, extract, io, contamination, cli


def test_paired_stat_empty():
    n, ngt, p, fold = stats.paired_stat([], [])
    assert n == 0 and ngt == 0
    assert math.isnan(p) and math.isnan(fold)


def test_paired_stat_all_tied_p_is_nan():
    a = np.array([1.0, 2.0, 3.0])
    n, ngt, p, fold = stats.paired_stat(a, a)
    assert n == 3 and ngt == 0
    assert math.isnan(p)            # no nonzero differences -> undefined
    assert fold == 1.0


def test_paired_stat_zero_denominator_is_finite():
    A = np.array([5.0, 3.0, 4.0])
    B = np.array([1.0, 0.0, 2.0])   # one zero denominator
    n, ngt, p, fold = stats.paired_stat(A, B)
    assert n == 3 and ngt == 3
    assert math.isfinite(fold)      # the inf ratio is dropped, not propagated


def test_paired_stat_direction_and_significance():
    A = np.array([10.0, 8.0, 9.0, 7.0, 11.0, 6.0])
    B = np.array([2.0, 1.5, 3.0, 1.0, 2.5, 1.2])
    n, ngt, p, fold = stats.paired_stat(A, B)
    assert ngt == 6 and fold > 1
    assert 0 < p <= 0.05


def test_extract_errors_on_empty_quantdir(tmp_path):
    cfg = {"groups": {"a": ["ENST1"], "b": ["ENST2"]}, "primary_comparison": ["a", "b"]}
    sm = tmp_path / "sm.csv"
    sm.write_text("donor,condition\nd1,control\n")
    with pytest.raises(FileNotFoundError):
        extract.extract(cfg, str(tmp_path / "nope"), str(sm), "C")


def test_primary_pair_requires_two_groups():
    with pytest.raises(ValueError):
        io.primary_pair({"groups": {"only": ["ENST1"]}, "primary_comparison": ["only"]})


def test_single_group_extract_still_writes(tmp_path):
    cfg = {"groups": {"only": ["ENST1"]}, "primary_comparison": ["only"]}
    sm = tmp_path / "sm.csv"; sm.write_text("donor,condition\nd1,control\n")
    d = tmp_path / "d1"; d.mkdir()
    (d / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nENST1.1\t100\t90\t5.0\t10\n")
    out = tmp_path / "pd.csv"
    n = extract.run(cfg, str(tmp_path), str(sm), "C", str(out))
    assert n == 1
    assert "only_TPM" in out.read_text()        # no pair fraction, no crash


def test_cli_kv_rejects_bad_argument():
    with pytest.raises(SystemExit):
        cli._kv(["GSE1=a.csv", "badarg"])          # missing '='
    assert cli._kv(["GSE1=a.csv"]) == {"GSE1": "a.csv"}


def test_contamination_missing_marker_columns_raises(tmp_path):
    p = tmp_path / "markers.csv"
    p.write_text("donor,SOMETHING\nd1,1.0\n")      # no tissue/contaminant columns
    with pytest.raises(ValueError):
        contamination.load_markers(str(p), ["TTR"], ["RBFOX3"])


def test_contamination_missing_target_column_raises(tmp_path):
    p = tmp_path / "target.csv"
    p.write_text("donor,other\nd1,1.0\n")
    with pytest.raises(ValueError):
        contamination.load_target(str(p), "long")   # no long_TPM column


def test_stats_empty_condition_raises(tmp_path):
    cfg = {"gene": "G", "groups": {"a": ["T1"], "b": ["T2"]}, "primary_comparison": ["a", "b"]}
    pd = tmp_path / "pd.csv"
    pd.write_text("cohort,donor,condition,a_TPM,b_TPM\nC,d1,control,5,1\n")
    with pytest.raises(ValueError):
        stats.run(cfg, "MISSING", {"C": str(pd)}, str(tmp_path / "res"))


def test_extract_errors_on_malformed_quantsf(tmp_path):
    cfg = {"groups": {"a": ["ENST1"], "b": ["ENST2"]}, "primary_comparison": ["a", "b"]}
    sm = tmp_path / "sm.csv"
    sm.write_text("donor,condition\nd1,control\n")
    d = tmp_path / "d1"
    d.mkdir()
    (d / "quant.sf").write_text("wrong\theader\n1\t2\n")   # no Name/TPM columns
    with pytest.raises(ValueError):
        extract.extract(cfg, str(tmp_path), str(sm), "C")
