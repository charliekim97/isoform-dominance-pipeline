"""End-to-end: synthetic data must reproduce the published LEPR result."""
import json
import re

from isoform_dominance import _selftest, cli


def test_reproduces_published_result(tmp_path):
    ok, msgs = _selftest.run(str(tmp_path))
    assert ok, "self-test failed:\n" + "\n".join(msgs)


def test_combined_pvalue(tmp_path):
    info = _selftest.generate(str(tmp_path))
    from isoform_dominance import extract, stats
    perdonor = {}
    for cohort, p in info.items():
        out = str(tmp_path / ("pd_%s.csv" % cohort))
        extract.run(_selftest.CONFIG, p["quantdir"], p["samplemap"], cohort, out)
        perdonor[cohort] = out
    res = stats.run(_selftest.CONFIG, "control", perdonor, str(tmp_path / "r"))
    cn, cgt, cp, _ = res["combined"]
    assert cn == 11 and cgt == 11
    assert abs(cp - 9.7656e-4) < 1e-4


def _reference_combination(tmp_path):
    """The three combinations, computed directly rather than through the self-test."""
    from isoform_dominance import extract, stats
    info = _selftest.generate(str(tmp_path))
    perdonor = {}
    for cohort, p in info.items():
        out = str(tmp_path / ("pd_%s.csv" % cohort))
        extract.run(_selftest.CONFIG, p["quantdir"], p["samplemap"], cohort, out)
        perdonor[cohort] = out
    return stats.run(_selftest.CONFIG, "control", perdonor, str(tmp_path / "r"),
                     n_boot=0)["combination"]


# --------------------------------------------------------------------------- #
# selftest --json
# --------------------------------------------------------------------------- #
def test_selftest_accepts_json_and_writes_only_json(capsys):
    """`cli.py` promises --json on every subcommand; selftest was the exception."""
    rc = cli.main(["selftest", "--json"])
    out = capsys.readouterr().out
    payload = json.loads(out)                 # the whole of stdout, so no stray text
    assert rc == 0
    assert payload["ok"] is True
    names = [c["name"] for c in payload["checks"]]
    assert names == ["GSE228458", "GSE137619", "COMBINED", "figure", "contamination_qc"]
    assert all(c["ok"] is True for c in payload["checks"])


def test_selftest_json_reports_a_failing_check_and_exits_one(monkeypatch, capsys):
    """`ok` must be computed, not constant: break one expectation and watch it flip."""
    monkeypatch.setitem(_selftest.EXPECT, "COMBINED", (11, 11, 0.5))
    rc = cli.main(["selftest", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert rc == 1
    assert payload["ok"] is False
    verdicts = {c["name"]: c["ok"] for c in payload["checks"]}
    assert verdicts.pop("COMBINED") is False
    assert all(verdicts.values())             # and only that one


def test_selftest_human_output_exit_code_is_unchanged(monkeypatch, capsys):
    assert cli.main(["selftest"]) == 0
    monkeypatch.setitem(_selftest.EXPECT, "COMBINED", (11, 11, 0.5))
    assert cli.main(["selftest"]) == 1
    assert "FAIL" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# all three cohort combinations reach the reader
# --------------------------------------------------------------------------- #
_COMBINATION_LINES = {"pooled": "donor-pooled", "stouffer": "Stouffer",
                      "stratified_signed_rank": "stratified signed-rank"}


def test_selftest_prints_all_three_combinations(tmp_path, capsys):
    """paper.md says the self-test reports the pooled, Stouffer and stratified
    combinations. It printed only the pooled one; the other two lived in a CSV in a
    temporary directory that was deleted on exit."""
    ref = _reference_combination(tmp_path)
    ps = [ref[k]["p"] for k in _COMBINATION_LINES]
    assert len({"%.4g" % p for p in ps}) == 3     # distinct, so no line can borrow another's

    assert cli.main(["selftest"]) == 0
    lines = capsys.readouterr().out.splitlines()
    for key, label in _COMBINATION_LINES.items():
        hits = [ln for ln in lines if re.match(r"\s*%s\s" % re.escape(label), ln)]
        assert len(hits) == 1, (label, lines)
        assert "P=%.4g" % ref[key]["p"] in hits[0], hits[0]


def test_selftest_json_carries_all_three_combinations(tmp_path, capsys):
    ref = _reference_combination(tmp_path)
    assert cli.main(["selftest", "--json"]) == 0
    got = json.loads(capsys.readouterr().out)["combinations"]
    assert set(got) == set(_COMBINATION_LINES)
    for key in _COMBINATION_LINES:
        assert got[key]["p"] == ref[key]["p"]
        assert got[key]["k"] == 2
