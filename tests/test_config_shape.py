"""A config whose classes or comparison have the wrong shape is refused, naming the place.

2.8.0 checked only that a config is a JSON object with a ``groups`` object.  On a small
synthetic config run with ``identifiability --sequences ... --no-gene-background``:

1. an empty class, ``{"A": [], "B": ["B1"]}``, exited 0 with a verdict
   (``weakly_identifiable``) -- a wrong answer, with nothing said;
2. a class that is a string, ``"A": "A1"``, was read one character at a time, as the
   transcript ids ``A`` and ``1``, which were then asked of Ensembl;
3. an id that is not a string, ``"A": [5]``, was an ``AttributeError`` traceback;
4. an empty id, ``"A": [""]``, was asked of Ensembl;
5. ``"primary_comparison": null`` and 6. ``"primary_comparison": 5`` were ``TypeError``
   tracebacks;
7. ``"primary_comparison": "AB"`` was read as ``["A", "B"]`` and exited 0.

``extract`` and ``stats`` read the same config through the same ``io.load_config`` and took
these as they came: tracebacks, or a table and exit 0.  Now the rule
(:func:`isoform_dominance.io.groups_shape_error`) is one line naming the config and the
place on the command line, exit 1, and a ``ValueError`` with the same sentence from
:func:`isoform_dominance.identifiability.analyze`.  The network is refused at once here: a
socket error would be retried with back-off, about 31 s a request.
"""
import json
import pathlib
import shutil
import urllib.request

import pytest

from isoform_dominance import cli, identifiability, io

ROOT = pathlib.Path(__file__).resolve().parents[1]
EXTRACT_ORDER = ROOT / "tests" / "data" / "extract_order"

SHARED = "ACGTTGCA" * 40
SEQS = {"A1": SHARED + "GATTACA" * 30, "B1": SHARED + "TACCATG" * 30}


def _cases(a, b, a_ids, b_ids):
    """The seven shapes, for a config whose classes are ``a`` and ``b``: ``{name: (groups,
    primary_comparison, the place the message must name)}``."""
    good = {a: a_ids, b: b_ids}
    return {
        "an empty class": ({a: [], b: b_ids}, [a, b], 'group "%s"' % a),
        "a class that is a string": ({a: a_ids[0], b: b_ids}, [a, b], 'group "%s"' % a),
        "an id that is not a string": ({a: [5], b: b_ids}, [a, b], 'group "%s"' % a),
        "an empty id": ({a: [""], b: b_ids}, [a, b], 'group "%s"' % a),
        "primary_comparison null": (good, None, "primary_comparison"),
        "primary_comparison a number": (good, 5, "primary_comparison"),
        "primary_comparison a string": (good, a + b, "primary_comparison"),
    }


CASES = _cases("A", "B", ["A1"], ["B1"])
REAL = _cases("short", "long", ["ENST00000000003"], ["ENST00000000001", "ENST00000000002"])


@pytest.fixture
def no_network(monkeypatch):
    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)


def _write(path, gene, groups, primary):
    path.write_text(json.dumps({"gene": gene, "groups": groups,
                                "primary_comparison": primary}))
    return str(path)


def _one_error_line(capsys, *needles):
    out, err = capsys.readouterr()
    assert out == ""
    lines = err.strip().splitlines()
    assert len(lines) == 1 and lines[0].startswith("error: "), err
    for n in needles:
        assert n in lines[0], (n, lines[0])
    return lines[0]


@pytest.mark.parametrize("case", sorted(CASES))
def test_identifiability_refuses_the_shape_in_one_line_and_asks_nothing(
        no_network, tmp_path, capsys, case):
    groups, primary, place = CASES[case]
    cfg = _write(tmp_path / "cfg.json", "G", groups, primary)
    sq = tmp_path / "seqs.json"
    sq.write_text(json.dumps(SEQS))
    assert cli.main(["identifiability", "--config", cfg, "--sequences", str(sq),
                     "--no-gene-background"]) == 1
    # the place names the class or the field; the config, the file
    _one_error_line(capsys, "config %s: " % cfg, place)


@pytest.mark.parametrize("command", ["extract", "stats"])
@pytest.mark.parametrize("case", sorted(REAL))
def test_extract_and_stats_refuse_the_same_shapes(tmp_path, capsys, command, case):
    # real inputs that run with the fixture's own config, so that what refuses is the
    # config's shape and nothing else
    groups, primary, place = REAL[case]
    work = tmp_path / "case"
    shutil.copytree(EXTRACT_ORDER, work)
    cfg = _write(work / "bad.json", "FAKE", groups, primary)
    argv = {"extract": ["--quantdir", str(work / "quant"), "--samplemap",
                        str(work / "samplemap.csv"), "--cohort", "C",
                        "--out", str(work / "perdonor.csv")],
            "stats": ["--perdonor", "C=%s" % (work / "perdonor_v230.csv"), "--n-boot", "0",
                      "--out", str(work / "fig")]}[command]
    good = [command, "--config", str(work / "config.json")] + argv
    assert cli.main(good) == 0
    capsys.readouterr()
    assert cli.main([command, "--config", cfg] + argv) == 1
    _one_error_line(capsys, "config %s: " % cfg, place)


@pytest.mark.parametrize("case", sorted(CASES))
def test_analyze_raises_a_value_error_with_the_same_sentence(no_network, tmp_path, case):
    groups, primary, place = CASES[case]
    with pytest.raises(ValueError) as raised:
        identifiability.analyze({"gene": "G", "groups": groups, "primary_comparison": primary},
                                sequences=SEQS, background_gene_transcripts=False)
    assert type(raised.value) is ValueError
    said = str(raised.value)
    assert said.startswith(place)
    cfg = _write(tmp_path / "cfg.json", "G", groups, primary)
    with pytest.raises(io.InputError) as loaded:
        io.load_config(cfg)
    assert str(loaded.value) == "config %s: %s" % (cfg, said)


def test_a_library_call_may_give_tuples(no_network):
    res = identifiability.analyze({"groups": {"A": ("A1",), "B": ("B1",)},
                                   "primary_comparison": ("A", "B")},
                                  sequences=SEQS, background_gene_transcripts=False)
    assert res["primary_comparison"] == ["A", "B"]
    assert res["groups"]["A"]["n_transcripts"] == 1


def test_a_comparison_may_be_left_out_and_qc_reads_no_groups(tmp_path):
    cfg = _write(tmp_path / "cfg.json", "G", {"A": ["A1"], "B": ["B1"]}, None)
    doc = json.loads(pathlib.Path(cfg).read_text())
    del doc["primary_comparison"]
    pathlib.Path(cfg).write_text(json.dumps(doc))
    assert io.load_config(cfg)["groups"] == {"A": ["A1"], "B": ["B1"]}
    # qc reads contamination_qc alone, and its config is not held to the rule
    bad = _write(tmp_path / "qc.json", "G", {"A": []}, "AB")
    assert io.load_config(bad, need_groups=False)["primary_comparison"] == "AB"


def _configs():
    """Every config in the repository: example/, the package's, the test fixtures' and the
    configs `annotate` wrote that the tests hold."""
    found = sorted(set(ROOT.glob("example/**/*config*.json"))
                   | set(ROOT.glob("src/isoform_dominance/data/*config*.json"))
                   | set(ROOT.glob("tests/data/**/config.json"))
                   | set(ROOT.glob("tests/data/v260_outputs/annotate_*.json")))
    assert len(found) >= 20, found
    return found


@pytest.mark.parametrize("path", _configs(), ids=lambda p: str(p.relative_to(ROOT)))
def test_every_config_in_the_repository_still_loads(path):
    cfg = io.load_config(str(path))
    assert io.groups_shape_error(cfg) is None
