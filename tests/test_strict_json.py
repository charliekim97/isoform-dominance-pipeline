"""Every JSON the package writes is standard JSON: no ``NaN``, ``Infinity`` or ``-Infinity``
(issue #20).

Python's ``json`` writes those tokens by default and reads them back, but RFC 8259 has no
such values and a standard parser (``jq``, JavaScript's ``JSON.parse``, R's ``jsonlite``)
refuses the whole document.  Through 2.6.0 the identifiability report of NTRK2 at release
116 carried nine ``Infinity``, and ``stats --json`` a ``NaN`` in every run.  A number that
is not finite is now written as ``null``, and the field beside it that says why --
``estimable``, ``finite_se``, ``defined``, ``fold_ci_method`` -- is checked here.

Reading stays lenient: a file an earlier version wrote is still read.
"""
import csv
import json
import pathlib
import random

import pytest

from isoform_dominance import _selftest, cli, io

DATA = pathlib.Path(__file__).parent / "data"
GTF = DATA / "gencode_mini" / "gencode.v50.mini.gtf.gz"
V260 = DATA / "v260_outputs"
K = 15


def _reject(token):
    raise AssertionError("non-standard JSON constant %s" % token)


def strict(text):
    """``json.loads`` as a standard parser reads: NaN, Infinity and -Infinity refused."""
    return json.loads(text, parse_constant=_reject)


def _rand(r, n):
    return "".join(r.choice("ACGT") for _ in range(n))


# --------------------------------------------------------------------------- #
# identifiability: a contrast estimable while its class totals are not
# --------------------------------------------------------------------------- #
def _not_estimable_case(tmp_path):
    """Two configured transcripts whose class totals leave the row space once the FASTA
    records that hold their junction windows are left out (--max-window-records 0), while
    their difference stays in it.  2.6.0 wrote 9 Infinity and 2 NaN for it, as many
    Infinity as NTRK2's report at release 116."""
    r = random.Random(10)
    p1 = _rand(r, 250)
    p2 = _rand(r, 2 * (250 - K + 1) + K - 1)
    seqs = {"ENST00000000001": p1, "ENST00000000002": p2}
    gene = {"ENST00000000011": p1 + p1 + p2}
    records = {"ENST90000000001": _rand(r, 40) + p1[-(K - 1):] + p1[:K - 1] + _rand(r, 40),
               "ENST90000000002": _rand(r, 40) + p1[-(K - 1):] + p2[:K - 1] + _rand(r, 40)}
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"gene": "GENEX", "groups": {"A": ["ENST00000000001"],
                                                           "B": ["ENST00000000002"]},
                               "primary_comparison": ["A", "B"]}))
    sq, bg, fa = tmp_path / "s.json", tmp_path / "b.json", tmp_path / "bg.fa"
    sq.write_text(json.dumps(seqs))
    bg.write_text(json.dumps(gene))
    fa.write_text("".join(">%s.1\n%s\n" % (t, s) for t, s in records.items()))
    return ["identifiability", "--config", str(cfg), "--sequences", str(sq),
            "--background-sequences", str(bg), "--background-fasta", str(fa),
            "--max-window-records", "0", "--k", str(K), "--window", str(K),
            "--min-log2fc", "0.5", "--json"]


def test_identifiability_json_is_standard_and_says_why_a_figure_is_null(tmp_path, capsys):
    rc = cli.main(_not_estimable_case(tmp_path))
    out = capsys.readouterr().out
    assert rc == cli.EXIT_EFFECT_NOT_RESOLVED
    rep = strict(out)
    for g in ("A", "B"):
        cls = rep["groups"][g]
        assert cls["estimable"] is False                       # says why ...
        assert cls["conditioning_factor"] is None              # ... these are null
        assert cls["gls_relative_se"] is None and cls["min_resolvable_log2fc"] is None
        assert cls["finite_se"] is False
    c = rep["contrast"]
    # estimable, and still no finite figure: `estimable` cannot say why, `finite_se` does
    assert c["estimable"] is True and c["finite_se"] is False
    assert c["conditioning_factor"] == pytest.approx(2 ** 0.5)
    assert c["gls_relative_se"] is None and c["min_resolvable_log2fc"] is None
    noise = rep["counting_noise"]
    assert noise["defined"] is False
    assert noise["log2_ratio_se"] is None and noise["min_resolvable_log2fc"] is None


def test_a_finite_figure_says_so(tmp_path, capsys):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"groups": {"A": ["A1"], "B": ["B1"]},
                               "primary_comparison": ["A", "B"]}))
    r = random.Random(3)
    shared = _rand(r, 1500)
    sq = tmp_path / "s.json"
    sq.write_text(json.dumps({"A1": shared + _rand(r, 800), "B1": shared + _rand(r, 1000)}))
    assert cli.main(["identifiability", "--config", str(cfg), "--sequences", str(sq),
                     "--json"]) == 0
    rep = strict(capsys.readouterr().out)
    assert rep["contrast"]["finite_se"] is True
    assert rep["counting_noise"]["defined"] is True
    assert all(rep["groups"][g]["finite_se"] for g in ("A", "B"))


def test_save_inputs_writes_standard_json(tmp_path, capsys):
    argv = _not_estimable_case(tmp_path)
    saved = tmp_path / "inputs.json"
    cli.main(argv + ["--save-inputs", str(saved)])
    capsys.readouterr()
    strict(saved.read_text())


# --------------------------------------------------------------------------- #
# stats, qc, extract, annotate, selftest
# --------------------------------------------------------------------------- #
def _perdonor(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cohort", "donor", "condition", "short_TPM", "long_TPM"])
        for donor, a, b in rows:
            w.writerow(["C", donor, "control", a, b])
    return str(path)


@pytest.mark.parametrize("n_boot", ["0", "200"])
def test_stats_json_is_standard(tmp_path, capsys, n_boot):
    """``combination.pooled.z`` is NaN in every run -- the pooled test is no z-combination
    -- and with --n-boot 0 every interval is NaN, which ``fold_ci_method`` None says."""
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"gene": "G", "groups": {"short": ["T1"], "long": ["T2"]},
                               "primary_comparison": ["short", "long"]}))
    a = _perdonor(tmp_path / "a.csv", [("d%d" % i, 5.0 + i, 1.0) for i in range(6)])
    assert cli.main(["stats", "--config", str(cfg), "--perdonor", "a=%s" % a,
                     "--n-boot", n_boot, "--out", str(tmp_path / "out"), "--json"]) == 0
    res = strict(capsys.readouterr().out)
    assert res["combination"]["pooled"]["z"] is None
    det = res["detail"][0]
    if n_boot == "0":
        assert det["fold_ci"] == [None, None] and det["fold_ci_method"] is None
    else:
        assert all(isinstance(x, float) for x in det["fold_ci"])
        assert det["fold_ci_method"] == "bootstrap"


@pytest.mark.filterwarnings("ignore::scipy.stats.ConstantInputWarning")
def test_qc_json_is_standard_when_a_correlation_is_undefined(tmp_path, capsys):
    """Constant marker scores: Spearman's rho is NaN."""
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(dict(_selftest.CONFIG)))
    markers = tmp_path / "m.csv"
    markers.write_text("donor,TTR,RBFOX3\n" + "".join("d%d,100,1\n" % i for i in range(5)))
    target = tmp_path / "t.csv"
    target.write_text("donor,long_TPM\n" + "".join("d%d,%d\n" % (i, i + 1) for i in range(5)))
    assert cli.main(["qc", "--config", str(cfg), "--markers", "C=%s" % markers,
                     "--target", "C=%s" % target, "--out", str(tmp_path / "qc"),
                     "--json"]) == 0
    rows = strict(capsys.readouterr().out)
    assert rows[0]["rho"] is None


def test_extract_json_and_its_files_are_standard(tmp_path, capsys):
    """extract is left as it is: what it writes is finite already, and stays so."""
    info = _selftest.generate(str(tmp_path / "work"))
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(_selftest.CONFIG))
    for cohort, p in info.items():
        out = tmp_path / ("pd_%s.csv" % cohort)
        assert cli.main(["extract", "--config", str(cfg), "--quantdir", p["quantdir"],
                         "--samplemap", p["samplemap"], "--cohort", cohort,
                         "--out", str(out), "--json"]) == 0
        strict(capsys.readouterr().out)
        strict(pathlib.Path(str(out) + ".index.json").read_text())
        with open(out) as f:
            for row in csv.DictReader(f):
                for key in ("short_TPM", "long_TPM", "short_fraction"):
                    assert float(row[key]) == float(row[key]) and abs(float(row[key])) < 1e300


def test_annotate_json_and_config_are_standard(tmp_path, capsys):
    out = tmp_path / "cfg.json"
    assert cli.main(["annotate", "--gene", "LEPR", "--gtf", str(GTF), "--out", str(out),
                     "--json"]) == 0
    assert strict(capsys.readouterr().out) == strict(out.read_text())


def test_selftest_json_is_standard(capsys):
    assert cli.main(["selftest", "--json"]) == 0
    strict(capsys.readouterr().out)


def test_a_writer_refuses_rather_than_writes_a_non_standard_token(tmp_path):
    """allow_nan=False is the guard: a non-finite number that reaches a writer unconverted
    is an error, not a file a standard parser refuses."""
    with pytest.raises(ValueError):
        io.dump_json({"x": float("nan")}, str(tmp_path / "x.json"))
    assert io.finite_json({"x": [float("inf"), -float("inf"), 1.5, (2.0, float("nan"))],
                           "y": True, "z": 3}) == {"x": [None, None, 1.5, [2.0, None]],
                                                   "y": True, "z": 3}


# --------------------------------------------------------------------------- #
# reading stays lenient
# --------------------------------------------------------------------------- #
def test_inputs_saved_by_2_6_0_are_read_and_rerun(tmp_path, capsys):
    """``identifiability_inputs_GENEX.json`` was written by the 2.6.0 release code
    (tests/data/v260_outputs/README.md); its report read |log2FC| 0.4808 for the contrast."""
    doc = json.loads((V260 / "identifiability_inputs_GENEX.json").read_text())
    assert doc["package_version"] == "2.6.0"
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"gene": "GENEX", "gene_id": "ENSG00000000009",
                               "groups": {"A": ["ENST00000000001"],
                                          "B": ["ENST00000000002"]},
                               "primary_comparison": ["A", "B"]}))
    assert cli.main(["identifiability", "--config", str(cfg), "--inputs",
                     str(V260 / "identifiability_inputs_GENEX.json"), "--min-log2fc", "0.5",
                     "--json"]) == 0
    rep = strict(capsys.readouterr().out)
    assert rep["verdict"] == "identifiable"
    assert rep["contrast"]["min_resolvable_log2fc"] == pytest.approx(0.48081641836346,
                                                                     rel=1e-9)


def test_a_file_holding_nan_or_infinity_is_still_read(tmp_path):
    """Hand-edited, or from a tool that writes Python's tokens: read as Python reads it."""
    p = tmp_path / "cfg.json"
    p.write_text('{"groups": {"A": ["T1"]}, "_note": NaN, "_w": Infinity}')
    cfg = io.load_config(str(p))
    assert cfg["_note"] != cfg["_note"] and cfg["_w"] == float("inf")
