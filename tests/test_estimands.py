"""Two estimands: the class fragment total, and the length-corrected class abundance.

Through 2.8.0 `identifiability` asked one question: whether a class's expected fragment count,
the sum of ``theta_t`` over its transcripts, is estimable.  `extract` and `stats` compare
summed TPM, whose class ratio is that of the length-corrected abundance, the sum of
``theta_t / l_t`` (``l_t`` the effective length).  The two are different functionals of one
system, and the same null space can leave one estimable and not the other.  A direction the
system cannot see that moves abundance between two transcripts of one class that differ in
length keeps the class's abundance and moves its fragment count by the difference in length,
so the fragment estimand is not estimable; the abundance one is when no direction the system
cannot see changes it.

A plain case: two alternative first exons ``X1``, ``X2`` crossed with two alternative last
exons ``Y1``, ``Y2`` across a shared middle ``M``, the classes by last exon.  Each class's
abundance is read off its own last exon, and is estimable; its fragment total moves along
``(+l_1, +l_2, -l_3, -l_4)`` when ``X1`` and ``X2`` differ in length, and is not.  With
every transcript of one length the two estimands are the same functional, and give the same
answer and the same figures.  With the classes taken across the cross instead,
``X1Y1 + X2Y2`` against ``X1Y2 + X2Y1``, neither is estimable: that direction moves both
classes' abundance.

`--estimand` picks the one that decides `verdict`, `reasons`, `effect_resolvable`, the figures
of `groups` and `contrast`, and the exit status: `fragment`, as 2.8.0 did, unless asked.  Both
are reported under `estimands`.
"""
import contextlib
import gzip
import io as _io
import json
import pathlib
import random

import pytest

from isoform_dominance import cli, identifiability as I

ROOT = pathlib.Path(__file__).resolve().parents[1]
SAVED = ROOT / "example" / "saved"
SAVED_RUNS = ["%s.gencode%d" % (g, n) for g in ("FLT1", "LEPR", "NTRK2", "NTRK3")
              for n in (44, 50)]
#: the fields of a class or of the contrast that depend on the estimand
ESTIMAND_FIELDS = ("estimable", "residual", "conditioning_factor", "rank", "label",
                   "gls_relative_se", "finite_se", "min_resolvable_log2fc", "beyond_linear",
                   "verdict", "reasons")


def _rand(r, n):
    return "".join(r.choice("ACGT") for _ in range(n))


def crossed(x1=300, x2=420, y1=350, y2=500, seed=7, exons=False):
    """Four transcripts: X1/X2 crossed with Y1/Y2 across a shared middle, the classes by last
    exon, ``Y1 = {T1, T4}`` and ``Y2 = {T2, T3}``; with ``exons``, the five exons too."""
    r = random.Random(seed)
    X1, X2, M, Y1, Y2 = (_rand(r, n) for n in (x1, x2, 700, y1, y2))
    seqs = {"T1": X1 + M + Y1, "T2": X2 + M + Y2, "T3": X1 + M + Y2, "T4": X2 + M + Y1}
    config = {"gene": "CROSS", "groups": {"Y1": ["T1", "T4"], "Y2": ["T2", "T3"]},
              "primary_comparison": ["Y1", "Y2"]}
    if exons:
        return config, seqs, {"X1": X1, "X2": X2, "M": M, "Y1": Y1, "Y2": Y2}
    return config, seqs


def _analyze(config, seqs, **kw):
    return I.analyze(config, sequences=seqs, background_gene_transcripts=False, **kw)


def _estimable(res, name):
    e = res["estimands"][name]
    return [e["groups"]["Y1"]["estimable"], e["groups"]["Y2"]["estimable"],
            e["contrast"]["estimable"]]


# --------------------------------------------------------------------------- #
# the algebra
# --------------------------------------------------------------------------- #
def test_crossed_classes_are_estimable_as_abundance_and_not_as_fragments():
    config, seqs = crossed()
    res = _analyze(config, seqs)
    assert _estimable(res, "fragment") == [False, False, False]
    assert _estimable(res, "abundance") == [True, True, True]
    # the abundance answer is exact, not a near miss of the threshold
    ab = res["estimands"]["abundance"]
    assert max(ab["groups"][g]["residual"] for g in ("Y1", "Y2")) < 1e-12
    assert ab["contrast"]["residual"] < 1e-12
    frag = res["estimands"]["fragment"]
    assert min(frag["groups"][g]["residual"] for g in ("Y1", "Y2")) > 1e-3
    # finite figures where the system resolves the estimand, none where it does not
    assert all(ab["groups"][g]["finite_se"] for g in ("Y1", "Y2"))
    assert ab["contrast"]["finite_se"]
    assert frag["contrast"]["min_resolvable_log2fc"] == float("inf")


def test_across_the_cross_neither_estimand_is_estimable():
    # the same sequences, the classes X1Y1 + X2Y2 against X1Y2 + X2Y1: crossed, of mixed
    # length, and estimable as neither, so crossing alone does not part the two
    config, seqs = crossed()
    config = dict(config, groups={"D1": ["T1", "T2"], "D2": ["T3", "T4"]},
                  primary_comparison=["D1", "D2"])
    res = _analyze(config, seqs)
    for name in I.ESTIMANDS:
        e = res["estimands"][name]
        assert [e["groups"]["D1"]["estimable"], e["groups"]["D2"]["estimable"],
                e["contrast"]["estimable"]] == [False, False, False], name


def test_the_effective_length_counts_every_window_position_dropped_or_not(tmp_path):
    """l_t is every window position of a transcript, as the denominator of its column is.
    With X1 and Y1 held by FASTA records and left out (--max-window-records 0), each class's
    abundance is still read off what is left -- Y1's off its junction with M -- and it is
    estimable only because l_t counts the dropped positions too."""
    config, seqs, exons = crossed(exons=True)
    fasta = tmp_path / "bg.fa"
    fasta.write_text(">ENST90000000001.1\n%s\n>ENST90000000002.1\n%s\n"
                     % (exons["X1"], exons["Y1"]))
    res = _analyze(config, seqs, background_fasta=str(fasta), max_window_records=0)
    dropped = res["background"]["windows_dropped"]["columns"]
    assert dropped == {"T1": 590, "T3": 270, "T4": 320}
    ab = res["estimands"]["abundance"]
    assert _estimable(res, "abundance") == [True, True, True]
    assert max(ab["groups"][g]["residual"] for g in ("Y1", "Y2")) < 1e-12
    assert ab["contrast"]["residual"] < 1e-12
    assert _estimable(res, "fragment") == [False, False, False]


def test_of_one_length_the_two_estimands_are_one_functional():
    # every transcript of one length: 1/l_t is one number, and the abundance functionals,
    # scaled to the fragment ones' norm, are the fragment ones
    config, seqs = crossed(x1=360, x2=360, y1=400, y2=400)
    res = _analyze(config, seqs)
    assert _estimable(res, "fragment") == [True, True, True]
    assert _estimable(res, "abundance") == [True, True, True]
    frag, ab = res["estimands"]["fragment"], res["estimands"]["abundance"]
    for key in ("min_resolvable_log2fc", "gls_relative_se", "conditioning_factor"):
        for g in ("Y1", "Y2"):
            assert ab["groups"][g][key] == pytest.approx(frag["groups"][g][key], rel=1e-9), key
        assert ab["contrast"][key] == pytest.approx(frag["contrast"][key], rel=1e-9), key
    assert ab["verdict"] == frag["verdict"]


def test_each_estimand_says_what_it_is():
    config, seqs = crossed()
    res = _analyze(config, seqs)
    ab, frag = res["estimands"]["abundance"], res["estimands"]["fragment"]
    assert ab["groups"]["Y1"]["label"] == "the Y1 class abundance"
    assert ab["contrast"]["label"] == "the class abundance contrast"
    # the fragment estimand's labels, and so its reasons, are 2.8.0's
    assert frag["groups"]["Y1"]["label"] == "the Y1 class total"
    assert frag["contrast"]["label"] == "the class contrast"
    assert "length-corrected abundance" in ab["definition"]
    assert "expected fragment count" in frag["definition"]
    # only the fragment estimand's precision figures were checked in simulation
    assert frag["precision_checked_in_simulation"] is True
    assert ab["precision_checked_in_simulation"] is False


def test_each_estimand_holds_only_its_own_fields():
    """`estimands` holds what depends on the estimand and nothing the report adds beside it,
    such as the contrast's effective-length fields, in the order the report gives them."""
    config, seqs = crossed()
    res = _analyze(config, seqs, min_log2fc=1.0)
    assert "class_mean_efflen" in res["contrast"]
    for name in I.ESTIMANDS:
        e = res["estimands"][name]
        assert list(e) == ["definition", "precision_checked_in_simulation", "groups",
                           "contrast", "effect_resolvable", "verdict", "reasons"]
        for g in ("Y1", "Y2"):
            assert list(e["groups"][g]) == list(ESTIMAND_FIELDS), (name, g)
        assert list(e["contrast"]) == list(ESTIMAND_FIELDS), name


# --------------------------------------------------------------------------- #
# which estimand decides
# --------------------------------------------------------------------------- #
def test_the_fragment_estimand_decides_unless_asked():
    config, seqs = crossed()
    res = _analyze(config, seqs, min_log2fc=5.0)
    assert res["decision_estimand"] == "fragment"
    assert res["verdict"] == res["estimands"]["fragment"]["verdict"] == "not_identifiable"
    assert res["effect_resolvable"] is False
    assert res["estimands"]["abundance"]["effect_resolvable"] is True
    res = _analyze(config, seqs, min_log2fc=5.0, estimand="abundance")
    assert res["decision_estimand"] == "abundance"
    assert res["effect_resolvable"] is True
    assert res["verdict"] == res["estimands"]["abundance"]["verdict"] != "not_identifiable"


@pytest.mark.parametrize("estimand", ["fragment", "abundance"])
def test_the_report_fields_are_the_deciding_estimands(estimand):
    config, seqs = crossed()
    res = _analyze(config, seqs, min_log2fc=1.0, estimand=estimand)
    chosen = res["estimands"][estimand]
    for g in ("Y1", "Y2"):
        assert {k: res["groups"][g][k] for k in ESTIMAND_FIELDS} == \
            {k: chosen["groups"][g][k] for k in ESTIMAND_FIELDS}
    assert {k: res["contrast"][k] for k in ESTIMAND_FIELDS} == \
        {k: chosen["contrast"][k] for k in ESTIMAND_FIELDS}
    assert (res["verdict"], res["reasons"], res["effect_resolvable"]) == \
        (chosen["verdict"], chosen["reasons"], chosen["effect_resolvable"])


def test_an_unknown_estimand_is_refused():
    config, seqs = crossed()
    with pytest.raises(ValueError, match="estimand must be one of fragment, abundance"):
        _analyze(config, seqs, estimand="tpm")


# --------------------------------------------------------------------------- #
# the command
# --------------------------------------------------------------------------- #
def _cli(tmp_path, *extra):
    config, seqs = crossed()
    cfg, sq = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps(config))
    sq.write_text(json.dumps(seqs))
    out, err = _io.StringIO(), _io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = cli.main(["identifiability", "--config", str(cfg), "--sequences", str(sq),
                         "--no-gene-background"] + list(extra))
    return code, out.getvalue(), err.getvalue()


def test_the_command_says_which_estimand_decides_and_reports_the_other(tmp_path):
    code, out, err = _cli(tmp_path, "--min-log2fc", "5")
    assert code == cli.EXIT_EFFECT_NOT_RESOLVED
    lines = out.splitlines()
    assert ("  estimand: fragment (class fragment totals) decides; abundance below the "
            "contrast (--estimand abundance)") in lines
    other = [x for x in lines if x.startswith("  abundance estimand: ")]
    assert len(other) == 1
    assert lines.index(other[0]) == 1 + next(i for i, x in enumerate(lines)
                                             if x.startswith("  contrast Y1 vs Y2: "))
    note = " ".join(err.split())
    assert ("Y1, Y2 and the class contrast are estimable as length-corrected class abundance "
            "-- whose class ratio is the ratio that summed TPM estimates, which `extract` and "
            "`stats` compare -- and not as fragment totals, which decide here: no direction "
            "this system cannot see changes any of them as abundance, while for each some "
            "direction changes it as a fragment total, by moving abundance between transcripts "
            "that differ in length.") in note
    assert "--estimand abundance decides by the abundance estimand" in note
    assert "have not yet been checked against a quantifier's TPM" in note
    assert "EFFECT SIZE: |log2FC| 5.000 NOT resolved" in out


def test_the_abundance_estimand_decides_when_asked(tmp_path):
    code, out, err = _cli(tmp_path, "--min-log2fc", "5", "--estimand", "abundance")
    assert code == cli.EXIT_OK
    lines = out.splitlines()
    assert ("  estimand: abundance (length-corrected class abundance, whose class ratio "
            "summed TPM estimates) decides; fragment below the contrast (--estimand fragment)"
            ) in lines
    assert any(x.startswith("  fragment estimand: [not_identifiable] Y1 min |log2FC| n/a")
               for x in lines)
    # what resolves the effect is the deciding estimand's
    assert ("  EFFECT SIZE: |log2FC| 5.000 resolved by both class abundances and the contrast "
            "at this design") in lines
    note = " ".join(err.split())
    assert "the abundance estimand decides (--estimand abundance)" in note
    assert "--estimand abundance decides by" not in note
    # where the two part is said whichever decides, and only the deciding one "decides here"
    assert ("Y1, Y2 and the class contrast are estimable as length-corrected class abundance "
            "-- whose class ratio is the ratio that summed TPM estimates, which `extract` and "
            "`stats` compare -- and not as fragment totals: no direction this system cannot "
            "see changes any of them as abundance") in note
    assert "which decide here" not in note
    # what its estimability rests on, and that its precision is not yet checked
    assert ("Its estimability is decided exactly in this window system, but where only it is "
            "estimable that answer rests on counting each transcript's length in window "
            "positions (identifiability.ESTIMANDS); its precision figures have not yet been "
            "checked against a quantifier's TPM in simulation") in note


def test_json_reports_both_and_names_the_deciding_one(tmp_path):
    code, out, _ = _cli(tmp_path, "--json")
    assert code == cli.EXIT_OK

    def refuse(token):
        raise AssertionError("non-standard JSON constant %s" % token)
    rep = json.loads(out, parse_constant=refuse)
    assert rep["decision_estimand"] == "fragment"
    assert set(rep["estimands"]) == {"fragment", "abundance"}
    # a figure that does not exist is null, with estimable and finite_se beside it
    c = rep["estimands"]["fragment"]["contrast"]
    assert c["min_resolvable_log2fc"] is None and c["finite_se"] is False
    assert c["estimable"] is False
    assert rep["estimands"]["abundance"]["precision_checked_in_simulation"] is False
    # after every key 2.8.0 wrote
    assert list(rep)[-2:] == ["decision_estimand", "estimands"]


def test_the_footnote_says_what_a_star_on_the_other_estimands_line_means(tmp_path):
    """SMN2 of the validation panel: as fragment totals, which decide, nothing is estimable and
    no figure is starred; as abundance its contrast's figure is past the linearisation limit,
    starred on the other estimand's line, and the footnote says what the star means."""
    entry = json.loads(gzip.open(ROOT / "validation" / "data" / "sweep116.json.gz",
                                 "rt", encoding="utf-8").read())["genes"]["SMN2"]
    config = entry["config"]
    seqs = {t.split(".")[0]: s for t, s in entry["sequences"].items()}
    named = {t.split(".")[0] for ids in config["groups"].values() for t in ids}
    paths = {}
    for key, part in (("config", config),
                      ("sequences", {t: s for t, s in seqs.items() if t in named}),
                      ("background-sequences", {t: s for t, s in seqs.items()
                                                if t not in named})):
        paths[key] = tmp_path / (key + ".json")
        paths[key].write_text(json.dumps(part))
    out, err = _io.StringIO(), _io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = cli.main(["identifiability"] + [x for key, path in paths.items()
                                                for x in ("--" + key, str(path))])
    assert code == cli.EXIT_OK
    lines = out.getvalue().splitlines()
    deciding = [x for x in lines if x.startswith(("  [", "  contrast "))]
    assert len(deciding) == 3 and not any("*" in x for x in deciding)
    other = [x for x in lines if x.startswith("  abundance estimand: ")]
    assert len(other) == 1 and "contrast min |log2FC| 0.867*" in other[0]
    assert ("  * first-order figure past the linearisation limit: read it as 'not resolvable "
            "at this design', not as a value.") in err.getvalue().splitlines()


def test_an_unknown_estimand_is_a_usage_error(capsys):
    with pytest.raises(SystemExit) as stop:
        cli.main(["identifiability", "--config", "x.json", "--estimand", "tpm"])
    assert stop.value.code == 2
    assert "invalid choice: 'tpm'" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# the saved examples: the default is 2.8.0's report, plus two keys
# --------------------------------------------------------------------------- #
def _saved(name, *extra):
    base = str(SAVED / name)
    out = _io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(_io.StringIO()):
        code = cli.main(["identifiability", "--config", base + ".config.json",
                         "--inputs", base + ".inputs.json", "--json"] + list(extra))
    return code, out.getvalue()


@pytest.mark.parametrize("name", SAVED_RUNS)
def test_asking_for_the_fragment_estimand_is_the_default(name):
    assert _saved(name) == _saved(name, "--estimand", "fragment")


@pytest.mark.parametrize("name", SAVED_RUNS)
def test_which_saved_examples_each_estimand_finds_estimable(name):
    # on every saved example the abundance estimand is estimable; the fragment one is not
    # on three (the converse can hold elsewhere: test_the_converse_is_said_too)
    rep = json.loads(_saved(name)[1])
    frag, ab = rep["estimands"]["fragment"], rep["estimands"]["abundance"]
    entries = [(frag["groups"][g], ab["groups"][g]) for g in rep["primary_comparison"]]
    entries.append((frag["contrast"], ab["contrast"]))
    assert all(a["estimable"] for _, a in entries)
    # NTRK2 at both releases and NTRK3 at GENCODE 50: crossed classes of mixed length
    assert all(f["estimable"] for f, _ in entries) == (
        name not in ("NTRK2.gencode44", "NTRK2.gencode50", "NTRK3.gencode50"))


def test_ntrk2_at_gencode_50_is_resolved_as_abundance():
    code, out = _saved("NTRK2.gencode50", "--min-log2fc", "0.5", "--estimand", "abundance")
    rep = json.loads(out)
    assert code == cli.EXIT_OK and rep["effect_resolvable"] is True
    assert rep["contrast"]["min_resolvable_log2fc"] == pytest.approx(0.1333, abs=5e-4)
    code, _ = _saved("NTRK2.gencode50", "--min-log2fc", "0.5")
    assert code == cli.EXIT_EFFECT_NOT_RESOLVED


def test_the_converse_is_said_too(tmp_path):
    """Where windows are dropped a column no longer sums to one, and the converse can hold:
    test_strict_json's case, whose contrast is estimable as fragment totals and not as
    abundance.  The NOTE says so, as it says the other way."""
    from test_strict_json import _not_estimable_case
    argv = [a for a in _not_estimable_case(tmp_path) if a != "--json"]
    out, err = _io.StringIO(), _io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        cli.main(argv + ["--json"])
    rep = json.loads(out.getvalue())
    assert rep["estimands"]["fragment"]["contrast"]["estimable"] is True
    assert rep["estimands"]["abundance"]["contrast"]["estimable"] is False
    out, err = _io.StringIO(), _io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        cli.main(argv)
    assert ("the class contrast is estimable as fragment totals and not as length-corrected "
            "class abundance: no direction this system cannot see changes it as a fragment "
            "total, while some change it as abundance, by moving abundance between transcripts "
            "that differ in length.") in " ".join(err.getvalue().split())
