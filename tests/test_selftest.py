"""End-to-end: synthetic data must reproduce the reference LEPR result."""
import json
import pathlib
import re
import shutil
from importlib import resources
from types import SimpleNamespace

import pytest

from isoform_dominance import _selftest, cli


def test_reproduces_the_reference_result(tmp_path):
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
    assert names == ["GSE228458", "GSE137619", "COMBINED", "figure", "contamination_qc",
                     "identifiability"]
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


# --------------------------------------------------------------------------- #
# identifiability: the study's classes on the saved GENCODE 44 inputs
# --------------------------------------------------------------------------- #
ROOT = pathlib.Path(__file__).resolve().parents[1]
PACKAGED = ROOT / "src" / "isoform_dominance" / "data"
SAVED = ROOT / "example" / "saved"


def _listing(folder):
    """The names of the files in ``folder``, sorted; a hidden one (``.DS_Store``) is none."""
    return sorted(p.name for p in folder.iterdir() if not p.name.startswith("."))


def test_a_hidden_file_is_not_in_the_listing(tmp_path):
    """A Finder's .DS_Store in the folder is no file of the package's."""
    copy = tmp_path / "data"
    shutil.copytree(PACKAGED, copy)
    (copy / ".DS_Store").write_bytes(b"\0\0\0\1Bud1")
    assert _listing(copy) == _listing(PACKAGED) == [
        "LEPR.gencode44.inputs.json", "LEPR.gencode44.study.config.json"]


def test_the_packaged_inputs_are_the_saved_ones_byte_for_byte():
    names = _listing(PACKAGED)
    assert names == ["LEPR.gencode44.inputs.json", "LEPR.gencode44.study.config.json"]
    for name in names:
        assert (PACKAGED / name).read_bytes() == (SAVED / name).read_bytes(), name
    # the study's classes: the ones the rest of the self-test quantifies
    study = json.loads((PACKAGED / "LEPR.gencode44.study.config.json").read_text())
    assert (study["groups"], study["primary_comparison"]) == (
        _selftest.CONFIG["groups"], _selftest.CONFIG["primary_comparison"])


def test_selftest_checks_the_study_classes_on_the_saved_inputs(capsys):
    assert cli.main(["selftest", "--json"]) == 0
    check = {c["name"]: c for c in json.loads(capsys.readouterr().out)["checks"]}[
        "identifiability"]
    assert check["ok"] is True
    assert check["verdict"] == "identifiable"
    assert abs(check["min_resolvable_log2fc"] - 0.1290) <= 0.0005
    assert check["min_resolvable_log2fc"] == pytest.approx(0.12896080608295227, rel=1e-9)
    assert check["conditioning_factor"] == pytest.approx(2.746955026236445, rel=1e-9)
    assert check["expect"] == {"verdict": "identifiable", "min_resolvable_log2fc": 0.129,
                               "tolerance": 0.0005}
    assert cli.main(["selftest"]) == 0
    line = [x for x in capsys.readouterr().out.splitlines() if "identifiability" in x]
    assert line == ["  [OK] identifiability: LEPR short vs long, saved GENCODE 44 inputs: "
                    "identifiable, contrast min |log2FC| 0.1290 (expect identifiable, "
                    "0.1290 +- 0.0005)"]


@pytest.mark.parametrize("change", [{"verdict": "weakly_identifiable"},
                                    {"min_resolvable_log2fc": 0.1300}])
def test_the_identifiability_check_can_fail(monkeypatch, capsys, change):
    monkeypatch.setattr(_selftest, "IDENTIFIABILITY", dict(_selftest.IDENTIFIABILITY, **change))
    assert cli.main(["selftest", "--json"]) == 1
    verdicts = {c["name"]: c["ok"] for c in json.loads(capsys.readouterr().out)["checks"]}
    assert verdicts.pop("identifiability") is False
    assert all(verdicts.values())


# --------------------------------------------------------------------------- #
# the package's own files, missing or broken: one failed check, not a lost report
# --------------------------------------------------------------------------- #
def _broken_package_data(monkeypatch, tmp_path, inputs):
    """Point the self-test's package data at a copy under ``tmp_path`` whose saved inputs
    are ``inputs`` (bytes), or absent when ``inputs`` is None."""
    data = tmp_path / "data"
    data.mkdir()
    name = _selftest.IDENTIFIABILITY["config"]
    (data / name).write_bytes((PACKAGED / name).read_bytes())
    if inputs is not None:
        (data / _selftest.IDENTIFIABILITY["inputs"]).write_bytes(inputs)
    monkeypatch.setattr(_selftest, "resources", SimpleNamespace(
        files=lambda package: tmp_path, as_file=resources.as_file))


_BROKEN = {"absent": (None, "No such file or directory"),
           "not JSON": (b'{"format": ', "is not valid JSON"),
           "not saved inputs": (b'{"format": "something else"}', "is not a saved-inputs file")}


@pytest.mark.parametrize("case", sorted(_BROKEN))
def test_broken_package_data_fails_one_check_and_keeps_the_json(monkeypatch, tmp_path,
                                                               capsys, case):
    inputs, cause = _BROKEN[case]
    _broken_package_data(monkeypatch, tmp_path, inputs)
    rc = cli.main(["selftest", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)          # the whole of stdout, so no stray text
    assert rc == 1
    assert payload["ok"] is False
    checks = {c["name"]: c for c in payload["checks"]}
    assert list(checks) == ["GSE228458", "GSE137619", "COMBINED", "figure",
                            "contamination_qc", "identifiability"]
    failed = checks.pop("identifiability")
    assert failed["ok"] is False
    assert failed["message"].startswith("[FAIL] identifiability: ")
    assert cause in failed["message"] and cause in failed["error"]
    assert (failed["verdict"], failed["min_resolvable_log2fc"]) == (None, None)
    assert all(c["ok"] is True for c in checks.values())
    assert set(payload["combinations"]) == set(_COMBINATION_LINES)
    assert payload["combinations"]["stouffer"]["p"] == pytest.approx(0.004419, abs=5e-7)


def test_broken_package_data_keeps_the_human_report(monkeypatch, tmp_path, capsys):
    _broken_package_data(monkeypatch, tmp_path, None)
    assert cli.main(["selftest"]) == 1
    lines = capsys.readouterr().out.splitlines()
    assert sum(ln.startswith("  [OK] ") for ln in lines) == 5
    failed = [ln for ln in lines if ln.startswith("  [FAIL] ")]
    assert len(failed) == 1 and failed[0].startswith("  [FAIL] identifiability: ")
    assert "No such file or directory" in failed[0]
    for label in _COMBINATION_LINES.values():
        assert any(re.match(r"\s*%s\s" % re.escape(label), ln) for ln in lines), label
    assert lines[-1] == "FAIL"
