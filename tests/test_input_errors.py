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


# --------------------------------------------------------------------------- #
# D1: a file that is not UTF-8
# --------------------------------------------------------------------------- #
LATIN1 = '{"gene": "G\xe9NE", "groups": {}}'.encode("latin-1")


@pytest.mark.parametrize("flag", ["--config", "--sequences", "--inputs"])
def test_a_file_that_is_not_utf8_is_one_line_naming_it(tmp_path, capsys, flag):
    bad = tmp_path / "bad.json"
    bad.write_bytes(LATIN1)
    seqs = tmp_path / "seqs.json"
    seqs.write_text(json.dumps({"ENST00000000001": "ACGT" * 100,
                                "ENST00000000002": "TTGA" * 100}))
    argv = {"--config": ["--config", str(bad), "--sequences", str(seqs)],
            "--sequences": ["--config", _cfg(tmp_path), "--sequences", str(bad)],
            "--inputs": ["--config", _cfg(tmp_path), "--inputs", str(bad)]}[flag]
    assert cli.main(["identifiability"] + argv) == 1
    _one_line(capsys, str(bad), "UTF-8")


# --------------------------------------------------------------------------- #
# D2: gzip is told by its first two bytes, not by the file name
# --------------------------------------------------------------------------- #
def _fasta_case(tmp_path, name, gz):
    import gzip as _gz
    shared = "ACGTTGCA" * 60
    seqs = {"ENST00000000001": shared + "A" * 3 + "GATTACA" * 50,
            "ENST00000000002": shared + "C" * 3 + "TACCATG" * 50}
    sq = tmp_path / "seqs.json"
    sq.write_text(json.dumps(seqs))
    text = ">ENSTBG0001 other\n%s\n" % ("GATTACA" * 50)
    fa = tmp_path / name
    fa.write_bytes(_gz.compress(text.encode()) if gz else text.encode())
    return ["identifiability", "--config", _cfg(tmp_path), "--sequences", str(sq),
            "--no-gene-background", "--background-fasta", str(fa), "--json"]


@pytest.mark.parametrize("name,gz", [("plain.fa.gz", False), ("bg", True), ("bg.fa.GZ", True),
                                     ("bg.fa", False), ("bg.fa.gz", True)])
def test_a_background_fasta_is_read_whatever_its_name_says(tmp_path, capsys, name, gz):
    assert cli.main(_fasta_case(tmp_path, name, gz)) == 0
    res = json.loads(capsys.readouterr().out)
    # the record shares the long class's tail, so it took unique k-mers from it
    plain = tmp_path / "ref"
    plain.mkdir()
    assert cli.main(_fasta_case(plain, "ref.fa", False)) == 0
    ref = json.loads(capsys.readouterr().out)
    assert res["groups"]["long"]["n_unique_kmers"] == ref["groups"]["long"]["n_unique_kmers"]


# --------------------------------------------------------------------------- #
# D3: qc with a config that has no contamination_qc, or a partial one
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("qc,needle", [
    (None, "contamination_qc"),
    ({"target_group": "long"}, "marker_panels"),
    ({"target_group": "long", "marker_panels": {"tissue": ["TTR"]}}, "contaminant"),
    ({"marker_panels": {"tissue": ["TTR"], "contaminant": ["PTPRC"]}}, "target_group"),
])
def test_qc_on_a_config_without_its_section_is_one_line(tmp_path, capsys, qc, needle):
    # what `annotate` writes has no contamination_qc: README step 4 runs qc on it
    cfg = _cfg(tmp_path, **({} if qc is None else {"contamination_qc": qc}))
    markers = tmp_path / "m.csv"
    markers.write_text("donor,TTR,PTPRC\nD1,5,1\nD2,6,2\nD3,7,1\n")
    target = tmp_path / "t.csv"
    target.write_text("cohort,donor,condition,long_TPM,short_TPM\n"
                      "C,D1,control,3,1\nC,D2,control,4,1\nC,D3,control,5,2\n")
    assert cli.main(["qc", "--config", cfg, "--markers", "C=" + str(markers),
                     "--target", "C=" + str(target), "--out", str(tmp_path / "qc")]) == 1
    _one_line(capsys, needle)


def test_qc_markers_without_a_panel_column_is_one_line(tmp_path, capsys):
    cfg = _cfg(tmp_path, contamination_qc={
        "target_group": "long", "marker_panels": {"tissue": ["TTR"], "contaminant": ["PTPRC"]}})
    markers = tmp_path / "m.csv"
    markers.write_text("donor,TTR\nD1,5\nD2,6\nD3,7\n")
    target = tmp_path / "t.csv"
    target.write_text("cohort,donor,condition,long_TPM\nC,D1,control,3\nC,D2,control,4\n"
                      "C,D3,control,5\n")
    assert cli.main(["qc", "--config", cfg, "--markers", "C=" + str(markers),
                     "--target", "C=" + str(target), "--out", str(tmp_path / "qc")]) == 1
    _one_line(capsys, str(markers))
