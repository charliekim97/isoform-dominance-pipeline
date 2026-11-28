"""Prospective G2 conditions, provenance, and failure rules without quantification."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "_gate_g2", Path(__file__).resolve().parents[1] / "validation" / "gate_g2.py")
g2 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(g2)


@pytest.fixture
def plan():
    controls = [f"C{i:02}" for i in range(41)]
    genes = sorted(g2.PRIMARY + controls + g2.REPORT_ONLY)
    predictions = {}
    for name in g2.CONDITIONS:
        predictions[name] = {
            "fragment": {g: {"predicted_se": 1.0 if g in controls else None,
                             "beyond_linear": False} for g in genes},
            "abundance": {g: {"predicted_se": None if g in g2.REPORT_ONLY else 1.0,
                              "beyond_linear": g == "CANX"} for g in genes},
        }
    return g2._plan(genes, [f"outside{i:04}" for i in range(3793)], predictions)


def _runs(plan, sha="frozen"):
    runs = {}
    for name, case in plan["conditions"].items():
        design = dict(g2.DESIGN, inputs_sha256=g2.INPUTS_SHA256,
                      genes_simulated=plan["genes"], prediction_background="reference_index",
                      competitor_transcripts=plan["competitor_transcripts"],
                      abundance_ratio=case["abundance_ratio"], seed_offset=case["seed_offset"],
                      replicate_seeds=case["seeds"], g2=g2._provenance(sha, name),
                      reference={"fasta_sha256": "fasta", "index": {
                          "salmon_version": "2.7.0", "k": 31, "seq_hash": "seq",
                          "name_hash": "name"}})
        rows = [{"gene": g, "n": 30, "bias": 0.1, "realised_sd": 0.5,
                 "predicted_se": 999.0} for g in plan["genes"]]
        runs[name] = {"design": design, "salmon": {"version": "2.7.0"},
                      "competitor_fragments": [100] * 30, "genes_tpm": rows}
    return runs


def _row(run, gene):
    return next(r for r in run["genes_tpm"] if r["gene"] == gene)


def test_seeds_are_distinct_and_flat_is_not_an_equal_class_target(plan):
    seeds = [seed for case in plan["conditions"].values() for seed in case["seeds"]]
    old = {7919 * rep + 13 for rep in range(30)}
    assert len(seeds) == len(set(seeds)) == 90
    assert old.isdisjoint(seeds)
    assert plan["conditions"]["flat"]["abundance_ratio"] is None


def test_plan_refuses_an_eligibility_change_and_an_incorrect_digest(plan, tmp_path):
    path = tmp_path / "plan.json"
    content = json.dumps(plan).encode()
    path.write_bytes(content)
    digest = hashlib.sha256(content).hexdigest()
    assert g2._read_frozen(path, digest) == plan
    with pytest.raises(ValueError, match="frozen digest"):
        g2._read_frozen(path, "wrong")
    plan["conditions"]["a4"]["predictions"]["abundance"]["CANX"]["predicted_se"] = None
    with pytest.raises(ValueError, match="eligibility"):
        g2._check_plan(plan)


def test_grader_uses_frozen_se_and_keeps_report_only_genes_out_of_the_gate(plan):
    runs = _runs(plan)
    _row(runs["flat"], "SMN2")["bias"] = 100.0
    result = g2.grade(plan, "frozen", runs)
    assert result["passes"]
    term = result["conditions"]["flat"]
    assert term["report_only"]["SMN2"]["bias"] == 100.0
    assert term["primary"]["genes"]["CANX"]["beyond_linear"]


@pytest.mark.parametrize("group,field,count", [
    ("primary", "bias", 1),
    ("primary", "realised_sd", 2),
    ("control", "bias", 5),
    ("control", "realised_sd", 10),
])
def test_each_count_rule_can_fail_on_equality_to_se(plan, group, field, count):
    runs = _runs(plan)
    for gene in plan["groups"][group][:count]:
        _row(runs["a4"], gene)[field] = 1.0
    result = g2.grade(plan, "frozen", runs)
    assert not result["passes"]
    assert not result["conditions"]["a4"][group]["passes"]
    assert result["conditions"]["b4"]["passes"]


def test_a_gross_control_error_stops_even_when_control_counts_pass(plan):
    runs = _runs(plan)
    _row(runs["b4"], "C00")["bias"] = -2.0
    group = g2.grade(plan, "frozen", runs)["conditions"]["b4"]["control"]
    assert group["bias_below_se"] == 40
    assert group["sd_below_se"] == 41
    assert group["gross"] == ["C00"]
    assert not group["passes"]


def test_unusable_rows_fail_with_the_original_denominator(plan):
    runs = _runs(plan)
    rows = runs["a4"]["genes_tpm"]
    rows.remove(_row(runs["a4"], "CANX"))
    _row(runs["a4"], "C00")["n"] = 2
    result = g2.grade(plan, "frozen", runs)["conditions"]["a4"]
    assert len(result["primary"]["genes"]) == 6
    assert not result["primary"]["passes"]
    assert not result["control"]["genes"]["C00"]["bias_below_se"]
    assert not result["control"]["genes"]["C00"]["sd_below_se"]


def test_grading_requires_the_frozen_provenance_and_all_conditions(plan):
    runs = _runs(plan)
    with pytest.raises(ValueError, match="each of flat"):
        g2.grade(plan, "frozen", {"flat": runs["flat"]})
    wrong = copy.deepcopy(runs)
    wrong["a4"]["design"]["replicate_seeds"][0] += 1
    with pytest.raises(ValueError, match="replicate_seeds"):
        g2.grade(plan, "frozen", wrong)
    wrong = copy.deepcopy(runs)
    wrong["a4"]["design"]["reference"]["index"]["name_hash"] = "another index"
    with pytest.raises(ValueError, match="different recorded index"):
        g2.grade(plan, "frozen", wrong)
    runs["flat"]["genes_tpm"].append(runs["flat"]["genes_tpm"][0])
    with pytest.raises(ValueError, match="duplicate TPM"):
        g2.grade(plan, "frozen", runs)


def test_condition_runner_uses_the_frozen_predictions_and_seed_offset(plan, tmp_path, monkeypatch):
    monkeypatch.setattr(g2, "_check_inputs", lambda *_: None)
    monkeypatch.setattr(g2.harness, "_reference_identity", lambda *_: {
        "fasta_sha256": "verified", "index": {"salmon_version": "2.7.0"},
    })
    monkeypatch.setattr(g2.subprocess, "run", lambda *_, **__: SimpleNamespace(stdout="salmon 2.7.0\n"))
    captured = {}

    def record(argv, *, frozen):
        captured.update(argv=argv, frozen=frozen)

    monkeypatch.setattr(g2.harness, "main", record)
    g2.run_condition(plan, "digest", "b4", tmp_path / "inputs", tmp_path / "reference",
                     tmp_path / "index", tmp_path / "work", tmp_path / "result", 4)
    argv = captured["argv"]
    assert argv[argv.index("--seed-offset") + 1] == "3000009"
    assert argv[argv.index("--abundance-ratio") + 1] == "0.25"
    assert argv[argv.index("--replicates") + 1] == "30"
    assert captured["frozen"]["predictions"] == plan["conditions"]["b4"]["predictions"]
    assert captured["frozen"]["provenance"] == g2._provenance("digest", "b4")


def test_condition_runner_rejects_result_paths_inside_replicate_work(plan, tmp_path):
    work = tmp_path / "work"
    with pytest.raises(ValueError, match="outside the replicate work directory"):
        g2.run_condition(plan, "digest", "flat", tmp_path / "inputs", tmp_path / "reference",
                         tmp_path / "index", work, work / "rep0.truth.json", 1)
    assert not work.exists()


def test_condition_runner_refuses_different_read_generation_code(plan, tmp_path):
    plan["execution"]["code_sha256"]["simulator"] = "different"
    with pytest.raises(ValueError, match="running code"):
        g2.run_condition(plan, "digest", "flat", tmp_path / "inputs", tmp_path / "reference",
                         tmp_path / "index", tmp_path / "work", tmp_path / "result", 1)
