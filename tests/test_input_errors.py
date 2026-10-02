"""A file or directory named on the command line that cannot be used is one line and exit 1.

Each case here was a traceback in 2.4.0: an ``OSError`` with no file name, which the CLI
re-raises because a network error looks the same, or a ``KeyError`` from deep inside.
"""
import json

import pytest

from isoform_dominance import cli

GROUPS = {"long": ["ENST00000000001"], "short": ["ENST00000000002"]}


def _one_line(capsys, *needles):
    err = capsys.readouterr().err
    assert "Traceback" not in err
    lines = err.strip().splitlines()
    assert len(lines) == 1, err
    for n in needles:
        assert n in lines[0], (n, lines[0])
    return lines[0]


def _cfg(tmp_path, **extra):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(dict({"gene": "G", "groups": GROUPS,
                                  "primary_comparison": ["short", "long"]}, **extra)))
    return str(p)


def _quant(tmp_path, donors=("D1", "D2")):
    qd = tmp_path / "quant"
    for d in donors:
        (qd / d).mkdir(parents=True)
        (qd / d / "quant.sf").write_text("Name\tTPM\nENST00000000001.1\t3\n"
                                         "ENST00000000002.1\t1\n")
    return str(qd)


# --------------------------------------------------------------------------- #
# extract
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("make", ["empty", "absent"])
def test_a_quantdir_with_no_quant_sf_is_one_line(tmp_path, capsys, make):
    qd = tmp_path / "quant"
    if make == "empty":
        qd.mkdir()
    sm = tmp_path / "sm.csv"
    sm.write_text("donor,condition\nD1,control\n")
    assert cli.main(["extract", "--config", _cfg(tmp_path), "--quantdir", str(qd),
                     "--samplemap", str(sm), "--cohort", "C",
                     "--out", str(tmp_path / "pd.csv")]) == 1
    _one_line(capsys, str(qd))


def test_a_sample_map_without_a_donor_column_is_one_line(tmp_path, capsys):
    sm = tmp_path / "sm.csv"
    sm.write_text("sample,condition\nD1,control\nD2,control\n")
    assert cli.main(["extract", "--config", _cfg(tmp_path), "--quantdir", _quant(tmp_path),
                     "--samplemap", str(sm), "--cohort", "C",
                     "--out", str(tmp_path / "pd.csv")]) == 1
    _one_line(capsys, str(sm), "donor")


# --------------------------------------------------------------------------- #
# stats
# --------------------------------------------------------------------------- #
def _perdonor(tmp_path, header, rows):
    p = tmp_path / "pd.csv"
    p.write_text(header + "\n" + "".join(r + "\n" for r in rows))
    return str(p)


def test_a_per_donor_table_without_a_class_column_is_one_line(tmp_path, capsys):
    pd = _perdonor(tmp_path, "cohort,donor,condition,long_TPM",
                   ["C,D%d,control,%d" % (i, i) for i in range(1, 6)])
    assert cli.main(["stats", "--config", _cfg(tmp_path), "--perdonor", "C=" + pd,
                     "--out", str(tmp_path / "fig")]) == 1
    _one_line(capsys, pd, "short_TPM")


def test_a_condition_no_donor_has_is_one_line(tmp_path, capsys):
    pd = _perdonor(tmp_path, "cohort,donor,condition,long_TPM,short_TPM",
                   ["C,D%d,disease,%d,1" % (i, i + 1) for i in range(1, 6)])
    assert cli.main(["stats", "--config", _cfg(tmp_path), "--perdonor", "C=" + pd,
                     "--condition", "control", "--out", str(tmp_path / "fig")]) == 1
    _one_line(capsys, pd, "control")
