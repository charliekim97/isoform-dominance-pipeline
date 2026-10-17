"""What `identifiability` says on a saved input does not move: the report 2.7.0 gave on every
saved example of `example/saved/`, with the same exit status.

The files compared against are in ``tests/data/v270_identifiability`` (its README says how they
were written).  Keys, strings, booleans, integers and nulls must be the same; a real number
within a relative 1e-9 (absolute 1e-12), because a report holds residuals of about 1e-16 whose
last digits move with the BLAS build.  The only change allowed is a key this release adds.
"""
import copy
import json
import math
import pathlib

import pytest

from isoform_dominance import cli

ROOT = pathlib.Path(__file__).resolve().parents[1]
SAVED = ROOT / "example" / "saved"
GOLDEN = ROOT / "tests" / "data" / "v270_identifiability"

#: keys of the report that 2.8.0 adds; everything else is 2.7.0's
ADDED = ("unique_read_counting_se",)

#: golden file -> (config, saved inputs) under example/saved
RUNS = {"%s.gencode%d" % (g, n): ("%s.gencode%d.config.json" % (g, n),
                                  "%s.gencode%d.inputs.json" % (g, n))
        for g in ("FLT1", "LEPR", "NTRK2", "NTRK3") for n in (44, 50)}
RUNS["LEPR.gencode44.study"] = ("LEPR.gencode44.study.config.json", "LEPR.gencode44.inputs.json")

#: the readable part: verdict, contrast estimable, contrast min_resolvable_log2fc to 4 places
SUMMARY = {
    "FLT1.gencode44": ("identifiable", True, 0.2717),
    "FLT1.gencode50": ("weakly_identifiable", True, 0.3270),
    "LEPR.gencode44": ("weakly_identifiable", True, 0.3183),
    "LEPR.gencode50": ("weakly_identifiable", True, 0.2342),
    "NTRK2.gencode44": ("not_identifiable", False, None),
    "NTRK2.gencode50": ("not_identifiable", False, None),
    "NTRK3.gencode44": ("identifiable", True, 0.2949),
    "NTRK3.gencode50": ("not_identifiable", False, None),
    "LEPR.gencode44.study": ("identifiable", True, 0.1290),
}

REL, ABS = 1e-9, 1e-12


def differences(got, want, where="report"):
    """Where two parsed reports differ, one line each; empty when they agree."""
    if isinstance(want, dict):
        if not isinstance(got, dict):
            return ["%s: %r is not an object" % (where, got)]
        out = ["%s: keys %s only in the new report" % (where, sorted(set(got) - set(want)))
               ] if set(got) - set(want) else []
        out += ["%s: keys %s missing" % (where, sorted(set(want) - set(got)))
                ] if set(want) - set(got) else []
        for k in sorted(set(got) & set(want)):
            out += differences(got[k], want[k], "%s.%s" % (where, k))
        return out
    if isinstance(want, list):
        if not isinstance(got, list) or len(got) != len(want):
            return ["%s: %r is not %r" % (where, got, want)]
        return [d for i, (g, w) in enumerate(zip(got, want, strict=True))
                for d in differences(g, w, "%s[%d]" % (where, i))]
    if isinstance(want, float) and type(got) is float:
        return [] if math.isclose(got, want, rel_tol=REL, abs_tol=ABS) else [
            "%s: %r is not %r" % (where, got, want)]
    # strings, booleans, integers and null: exactly, and of the same type
    return [] if type(got) is type(want) and got == want else [
        "%s: %r is not %r" % (where, got, want)]


def _report(name, capsys):
    config, inputs = RUNS[name]
    code = cli.main(["identifiability", "--config", str(SAVED / config),
                     "--inputs", str(SAVED / inputs), "--json"])
    return code, json.loads(capsys.readouterr().out)


@pytest.mark.parametrize("name", sorted(RUNS))
def test_the_report_on_a_saved_input_is_the_one_2_7_0_gave(name, capsys):
    code, got = _report(name, capsys)
    want = json.loads((GOLDEN / (name + ".json")).read_text())
    assert code == 0
    for key in ADDED:
        got.pop(key, None)
    assert differences(got, want) == []
    verdict, estimable, fc = SUMMARY[name]
    c = got["contrast"]
    assert (got["verdict"], c["estimable"]) == (verdict, estimable)
    assert (c["min_resolvable_log2fc"] if fc is None
            else round(c["min_resolvable_log2fc"], 4)) == fc


def test_the_comparison_can_fail():
    """Each kind of change the comparison is there to catch makes it report one."""
    want = json.loads((GOLDEN / "LEPR.gencode50.json").read_text())
    assert differences(copy.deepcopy(want), want) == []
    changes = [
        lambda r: r["contrast"].__setitem__("min_resolvable_log2fc",
                                            r["contrast"]["min_resolvable_log2fc"] * (1 + 1e-7)),
        lambda r: r["contrast"].__setitem__("estimable", False),
        lambda r: r.__setitem__("verdict", "identifiable"),
        lambda r: r["groups"]["iso_896aa"].__setitem__("n_unique_kmers", 357),
        lambda r: r["groups"]["iso_896aa"].__setitem__("n_unique_kmers", 356.0),
        lambda r: r.__setitem__("new_key", 1),
        lambda r: r.pop("counting_noise"),
        lambda r: r["background"]["gene_transcripts"].pop(),
        lambda r: r["contrast"].__setitem__("gls_relative_se", None),
    ]
    for change in changes:
        got = copy.deepcopy(want)
        change(got)
        assert differences(got, want), change
