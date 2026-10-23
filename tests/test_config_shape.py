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
7. ``"primary_comparison": "AB"`` was read as ``["A", "B"]`` and exited 0;
8. no class at all, ``"groups": {}``, was a ``config error:`` line from ``identifiability``;
   with a comparison, a ``KeyError`` after one header line from ``extract`` and a run on the
   per-donor table from ``stats``, and without one a table of no class (exit 0) and a
   traceback;
9. an id with white space, ``"A": ["A1 "]``, was asked of Ensembl, and ``extract`` found it in
   no ``quant.sf``; with a version, ``"A1.1 "``, all three read it up to the dot, as if clean;
10. a comparison entry that is not a string, ``["A", 5]``, was a ``config error:`` line
    naming the missing group ``5`` from ``identifiability``;
11. a class given as an object, ``"A": {"A1": 1}``, was read as its keys by all three.

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

import numpy as np
import pytest

from isoform_dominance import cli, extract, identifiability, io

ROOT = pathlib.Path(__file__).resolve().parents[1]
EXTRACT_ORDER = ROOT / "tests" / "data" / "extract_order"

SHARED = "ACGTTGCA" * 40
SEQS = {"A1": SHARED + "GATTACA" * 30, "B1": SHARED + "TACCATG" * 30}


def _cases(a, b, a_ids, b_ids):
    """The twelve shapes, for a config whose classes are ``a`` and ``b``: ``{name: (groups,
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
        "a comparison entry not a string": (good, [a, 5], "primary_comparison lists 5"),
        "no class": ({}, [a, b], "groups is empty"),
        "an id with white space": ({a: [a_ids[0] + " "], b: b_ids}, [a, b], 'group "%s"' % a),
        "a versioned id with white space": ({a: [a_ids[0] + ".1 "], b: b_ids}, [a, b],
                                            'group "%s"' % a),
        "a class given as an object": ({a: dict.fromkeys(a_ids, 1), b: b_ids}, [a, b],
                                       'group "%s" is {' % a),
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


SEQS3 = dict(SEQS, A2=SHARED + "GGGCCCAT" * 25)


def _figures(res):
    # a set is iterated in Python's order, so the sums agree to the last digit or so
    return (res["verdict"], list(res["primary_comparison"]),
            {g: (r["n_transcripts"], r["estimable"]) for g, r in res["groups"].items()},
            pytest.approx(res["contrast"]["min_resolvable_log2fc"], rel=1e-12),
            pytest.approx(res["contrast"]["gls_relative_se"], rel=1e-12))


@pytest.mark.parametrize("groups,primary", [
    ({"A": {"A1", "A2"}, "B": frozenset(["B1"])}, None),
    ({"A": np.array(["A1", "A2"]), "B": np.array(["B1"])}, np.array(["A", "B"])),
    ({"A": ["A1", "A2"], "B": ["B1"]}, ("A", "B")),
], ids=["sets", "arrays", "a tuple"])
def test_a_library_call_may_give_any_collection_2_8_0_took(no_network, groups, primary):
    # 2.8.0 read a class as whatever it could iterate, and a comparison as whatever it could
    # index; the shape rule refuses strings, mappings and empty classes, not these
    as_lists = identifiability.analyze({"groups": {"A": ["A1", "A2"], "B": ["B1"]},
                                        "primary_comparison": ["A", "B"]},
                                       sequences=SEQS3, background_gene_transcripts=False)
    cfg = {"groups": groups}
    if primary is not None:
        cfg["primary_comparison"] = primary
    res = identifiability.analyze(cfg, sequences=SEQS3, background_gene_transcripts=False)
    assert _figures(res) == _figures(as_lists)
    assert io.groups_shape_error(cfg) is None


def test_a_set_has_no_order_to_compare_in_and_a_mapping_is_no_class(no_network):
    with pytest.raises(ValueError, match="primary_comparison is .*, a set: it has no first"):
        identifiability.analyze({"groups": {"A": ["A1"], "B": ["B1"]},
                                 "primary_comparison": {"A", "B"}},
                                sequences=SEQS, background_gene_transcripts=False)
    with pytest.raises(ValueError, match='group "A" is .*, not a list of transcript ids'):
        identifiability.analyze({"groups": {"A": {"A1": 1}, "B": ["B1"]}},
                                sequences=SEQS, background_gene_transcripts=False)


def test_analyze_without_groups_names_them():
    # 2.8.0: KeyError: 'groups'
    with pytest.raises(ValueError, match="^groups is null, not an object"):
        identifiability.analyze({}, sequences=SEQS, background_gene_transcripts=False)


# --------------------------------------------------------------------------------------- #
# How many groups the comparison names, a group the config lacks and a group named twice are
# checked where the comparison is used.  `stats` uses it first: through 2.8.0 a comparison
# naming one group, or one group twice, was a ValueError traceback, and one naming a group
# the config lacks was refused by the per-donor table, for a column the config never named.
# `extract` writes the comparison's fraction column, and for a group the config lacks it
# raised KeyError.
@pytest.mark.parametrize("primary,said", [
    (["short"], "primary_comparison must name two isoform groups"),
    (["short", "medium"], "primary_comparison names group(s) not in config['groups']: "
                          "['medium']"),
    (["short", "short"], "primary_comparison compares group 'short' with itself"),
    # 2.8.0 compared the first two, and ran
    (["short", "long", "short"], "primary_comparison must name two isoform groups for a "
                                 "paired comparison; got ['short', 'long', 'short'], 3 groups"),
    (["short", "long", "medium"], "primary_comparison must name two isoform groups"),
])
def test_stats_refuses_a_comparison_it_cannot_make_in_one_line(tmp_path, capsys, primary,
                                                                said):
    work = tmp_path / "case"
    shutil.copytree(EXTRACT_ORDER, work)
    cfg = _write(work / "bad.json", "FAKE",
                 {"long": ["ENST00000000001", "ENST00000000002"],
                  "short": ["ENST00000000003"]}, primary)
    argv = ["stats", "--config", cfg, "--perdonor", "C=%s" % (work / "perdonor_v230.csv"),
            "--n-boot", "0", "--out", str(work / "fig")]
    assert cli.main(argv) == 1
    _one_error_line(capsys, "error: config %s: %s" % (cfg, said))
    assert not (work / "fig_stats.csv").exists()


@pytest.mark.parametrize("primary,says", [
    # through 2.8.0 a KeyError, after one header line
    (["short", "medium"], 'names "medium", which the config\'s groups do not have'),
    # through 2.8.0 no fraction column, and nothing said
    (["short", "long", "short"], "names 3 group(s), not two"),
    (["short"], "names 1 group(s), not two"),
])
def test_extract_says_why_it_writes_no_fraction(tmp_path, capsys, primary, says):
    work = tmp_path / "case"
    shutil.copytree(EXTRACT_ORDER, work)
    cfg = _write(work / "bad.json", "FAKE",
                 {"long": ["ENST00000000001", "ENST00000000002"],
                  "short": ["ENST00000000003"]}, primary)
    out = work / "perdonor.csv"
    assert cli.main(["extract", "--config", cfg, "--quantdir", str(work / "quant"),
                     "--samplemap", str(work / "samplemap.csv"), "--cohort", "C",
                     "--out", str(out)]) == 0
    assert out.read_text().splitlines()[0] == "cohort,donor,condition,long_TPM,short_TPM"
    warned = [x for x in capsys.readouterr().err.splitlines() if "primary_comparison" in x]
    assert warned == ["  WARNING: primary_comparison %s, so no fraction column is written; "
                      "`identifiability` and `stats` refuse such a config" % says]


def test_a_set_comparison_is_refused_where_it_is_read_in_order(tmp_path):
    # load_config and analyze refuse it by its shape; primary_pair (stats.run) and extract,
    # called from Python, read the comparison too, and list() of a set takes Python's order
    # of the moment: a fraction column and a stats direction that change with the hash seed
    cfg = {"groups": {"long": ["ENST00000000001", "ENST00000000002"],
                      "short": ["ENST00000000003"]}, "primary_comparison": {"short", "long"}}
    for call in (lambda: io.primary_pair(cfg),
                 lambda: extract.write_perdonor(cfg, [], "C", str(tmp_path / "pd.csv")),
                 lambda: extract.run(cfg, str(EXTRACT_ORDER / "quant"),
                                     str(EXTRACT_ORDER / "samplemap.csv"), "C",
                                     str(tmp_path / "pd2.csv"))):
        with pytest.raises(ValueError, match="primary_comparison is .*, a set: it has no first"):
            call()
    assert not (tmp_path / "pd.csv").exists() and not (tmp_path / "pd2.csv").exists()


def test_analyze_names_a_missing_group_once(no_network):
    with pytest.raises(ValueError, match=r"not in config\['groups'\]: \['C'\]$"):
        identifiability.analyze({"groups": {"A": ["A1"], "B": ["B1"]},
                                 "primary_comparison": ["C", "C"]},
                                sequences=SEQS, background_gene_transcripts=False)


def test_extract_writes_no_fraction_for_a_group_the_config_lacks(tmp_path, capsys):
    work = tmp_path / "case"
    shutil.copytree(EXTRACT_ORDER, work)
    cfg = _write(work / "bad.json", "FAKE",
                 {"long": ["ENST00000000001", "ENST00000000002"],
                  "short": ["ENST00000000003"]}, ["short", "medium"])
    out = work / "perdonor.csv"
    assert cli.main(["extract", "--config", cfg, "--quantdir", str(work / "quant"),
                     "--samplemap", str(work / "samplemap.csv"), "--cohort", "C",
                     "--out", str(out)]) == 0
    assert out.read_text().splitlines()[0] == "cohort,donor,condition,long_TPM,short_TPM"
    err = capsys.readouterr().err
    assert ('WARNING: primary_comparison names "medium", which the config\'s groups do not '
            "have, so no fraction column is written; `identifiability` and `stats` refuse "
            "such a config") in err


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
