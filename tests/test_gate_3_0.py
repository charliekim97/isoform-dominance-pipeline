"""The TPM truth of the harness and the grader of the 3.0 gate (``validation/README.md``).

``run_simulation.tpm_truth`` is the gate's truth: per class, the fragments written for each
transcript over the effective length the quantifier gives it, a transcript the index dropped
as identical to another taking the kept one's.  ``gate_3_0.py`` applies the gate's steps 3
to 6 to a run's TPM rows against the fixed predictions; its numbers are held here to the
gate's text, and every way it can grade is tried on rows made from those predictions.
"""
import copy
import hashlib
import importlib.util
import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"


def _module(name):
    spec = importlib.util.spec_from_file_location("_" + name, VALIDATION / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


harness = _module("run_simulation")
gate = _module("gate_3_0")
PREDICTIONS = json.loads((VALIDATION / "gate_3_0_predictions.json").read_text(encoding="utf-8"))
STORED = json.loads((VALIDATION / "data" / "sim" / "c1_d1.json").read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# the TPM truth
# --------------------------------------------------------------------------- #
GROUPS = {"A": ["a1", "a2", "a3"], "B": ["b1"]}
SEQS = {"a1": "ACGT" * 50, "a2": "ACGT" * 50, "a3": "TTGCA" * 60, "b1": "GGCAT" * 80}


def test_the_truth_is_the_fragments_over_the_quantifiers_effective_length():
    efflen = {"a1": 20.0, "a2": 20.0, "a3": 110.0, "b1": 230.0}
    written = {"a1": 4, "a2": 6, "a3": 11, "b1": 46}
    truth, took = harness.tpm_truth(GROUPS, ["A", "B"], SEQS, written, efflen)
    assert truth == pytest.approx({"A": 4 / 20 + 6 / 20 + 11 / 110, "B": 46 / 230})
    assert took == {}


def test_a_transcript_the_index_dropped_takes_the_kept_ones_length():
    # a2 is identical to a1 and the index kept a1: a2's fragments stay, over a1's length
    efflen = {"a1": 20.0, "a3": 110.0, "b1": 230.0}
    written = {"a1": 4, "a2": 6, "a3": 11, "b1": 46}
    truth, took = harness.tpm_truth(GROUPS, ["A", "B"], SEQS, written, efflen)
    assert truth["A"] == pytest.approx(10 / 20 + 11 / 110)
    assert took == {"a2": "a1"}
    # the index's duplicate_clusters.tsv names the kept one, even outside the gene
    efflen = {"x9": 25.0, "a3": 110.0, "b1": 230.0}
    seqs = dict(SEQS, a1="CCCC" * 50)            # no identical copy within the gene
    truth, took = harness.tpm_truth(GROUPS, ["A", "B"], seqs, written, efflen,
                                    kept={"a1": "x9", "a2": "x9"})
    assert truth["A"] == pytest.approx(4 / 25 + 6 / 25 + 11 / 110)
    assert took == {"a1": "x9", "a2": "x9"}


def test_a_truth_that_would_be_undefined_is_refused():
    written = {"a1": 4, "a2": 6, "a3": 11, "b1": 46}
    with pytest.raises(RuntimeError, match="lists no a3 and no transcript of identical"):
        harness.tpm_truth(GROUPS, ["A", "B"], SEQS, written,
                          {"a1": 20.0, "a2": 20.0, "b1": 230.0})
    with pytest.raises(RuntimeError, match="effective length the quantifier gives b1 is 0.0"):
        harness.tpm_truth(GROUPS, ["A", "B"], SEQS, written,
                          {"a1": 20.0, "a2": 20.0, "a3": 110.0, "b1": 0.0})


def test_the_kept_copies_are_read_from_the_index(tmp_path):
    assert harness.kept_copies(str(tmp_path)) == {}
    (tmp_path / "duplicate_clusters.tsv").write_text(
        "RetainedTxp\tDuplicateTxp\nENST01.2\tENST02.1\nENST01.2\tENST03.4\n")
    assert harness.kept_copies(str(tmp_path)) == {"ENST02": "ENST01", "ENST03": "ENST01"}


# --------------------------------------------------------------------------- #
# the grader: its numbers are the gate's
# --------------------------------------------------------------------------- #
def _text():
    return " ".join((VALIDATION / "README.md").read_text(encoding="utf-8").split())


def test_the_graders_numbers_are_the_gates():
    text = _text()
    assert ("every one of the 7 has a |bias| below its predicted SE, at least %d of the 7 a "
            "replicate SD below it" % gate.PRIMARY_SD_AT_LEAST) in text
    assert ("Two of the 7 already failed with `NumReads`: %s (|bias|" % gate.STOP_GENES[0]) in text
    assert "and %s (|bias|" % gate.STOP_GENES[1] in text
    assert ("a |bias| below the predicted SE for at least %d, and an SD below it for at least %d"
            % (gate.CONTROL_BIAS_AT_LEAST, gate.CONTROL_SD_AT_LEAST)) in text
    assert "Salmon %s against the stored index" % gate.SALMON_VERSION in text
    assert len(PREDICTIONS["primary"]) == 7 and len(PREDICTIONS["control"]) == 39
    assert set(gate.STOP_GENES) <= set(PREDICTIONS["primary"])
    stored = json.loads((VALIDATION / "data" / "sim" / "c1_d1.json").read_text())["design"]
    assert {k: v for k, v in gate.DESIGN.items() if k != "quantifier"} == stored


def test_the_predictions_are_checked_against_the_gates_sha256(tmp_path):
    assert gate.frozen_predictions() == PREDICTIONS
    changed = tmp_path / "p.json"
    changed.write_bytes((VALIDATION / "gate_3_0_predictions.json").read_bytes() + b" ")
    with pytest.raises(gate.NotGradeable, match="not the .* the gate names"):
        gate.frozen_predictions(str(changed))


# --------------------------------------------------------------------------- #
# the grader: every way it can grade
# --------------------------------------------------------------------------- #
def _run(scale=0.5, sd_scale=0.5, overrides=None, drop=()):
    """A run whose every row has a |bias| and an SD of ``scale`` and ``sd_scale`` times its
    predicted SE; ``overrides`` sets a gene's (bias, SD) as multiples of its SE."""
    overrides = overrides or {}
    rows = []
    for part in ("primary", "control"):
        for gene, se in PREDICTIONS[part].items():
            if gene in drop:
                continue
            b, s = overrides.get(gene, (scale, sd_scale))
            rows.append({"gene": gene, "n": 30, "bias": -b * se, "realised_sd": s * se,
                         "rmse": 0.0})
    design = dict(gate.DESIGN, genes_simulated=sorted(
        set(PREDICTIONS["primary"]) | set(PREDICTIONS["control"]) | {"EPOR", "PFKL"}))
    return {"design": design, "salmon": {"version": gate.SALMON_VERSION},
            "genes": copy.deepcopy(STORED["genes"]), "genes_tpm": rows}


def test_a_run_inside_every_bound_passes():
    r = gate.grade(_run(), PREDICTIONS)
    assert r["primary_passes"] and r["control"]["passes"] and r["passes"]
    assert r["stopped"] == []


def test_one_sd_of_the_primary_set_may_miss_but_not_smn2s_or_vegfas():
    assert gate.grade(_run(overrides={"CANX": (0.5, 1.0)}), PREDICTIONS)["passes"]
    for gene in gate.STOP_GENES:
        r = gate.grade(_run(overrides={gene: (0.5, 1.0)}), PREDICTIONS)
        assert r["stopped"] == [gene] and not r["primary_passes"] and not r["passes"]
        assert r["control"] is None                 # graded only when the primary set passes
    r = gate.grade(_run(overrides={"CANX": (0.5, 1.2), "GNAS": (0.5, 1.0)}), PREDICTIONS)
    assert r["primary_sd_ok"] == 5 and not r["passes"]


def test_a_bias_at_the_predicted_se_fails_the_primary_set():
    r = gate.grade(_run(overrides={"NTRK3": (1.0, 0.5)}), PREDICTIONS)
    assert not r["primary_passes"] and r["stopped"] == []
    r = gate.grade(_run(overrides={"VEGFA": (1.0, 0.5)}), PREDICTIONS)
    assert r["stopped"] == ["VEGFA"]


def test_a_gene_without_a_row_fails_both_criteria():
    r = gate.grade(_run(drop=("SMN2",)), PREDICTIONS)
    assert r["stopped"] == ["SMN2"] and not r["passes"]
    r = gate.grade(_run(drop=("ALDOA",)), PREDICTIONS)
    assert r["primary_passes"]
    assert (r["control"]["bias_ok"], r["control"]["sd_ok"]) == (38, 38) and r["passes"]


def test_the_control_needs_38_biases_and_31_sds():
    control = list(PREDICTIONS["control"])
    two_bias = dict.fromkeys(control[:2], (1.5, 0.5))
    r = gate.grade(_run(overrides=two_bias), PREDICTIONS)
    assert r["control"]["bias_ok"] == 37 and not r["passes"]
    eight_sd = dict.fromkeys(control[:8], (0.5, 1.5))
    r = gate.grade(_run(overrides=eight_sd), PREDICTIONS)
    assert r["control"]["sd_ok"] == 31 and r["passes"]
    nine_sd = dict.fromkeys(control[:9], (0.5, 1.5))
    assert not gate.grade(_run(overrides=nine_sd), PREDICTIONS)["passes"]


def test_a_run_the_gate_does_not_name_is_not_graded():
    good = _run()
    gate.check_run(good, PREDICTIONS)
    for change, why in (
            (lambda r: r["design"].__setitem__("bias_5p", 0.3), "bias_5p 0.3"),
            (lambda r: r["design"].__setitem__("replicates", 29), "replicates 29"),
            (lambda r: r["design"].__setitem__("salmon_flags", "--seqBias"), "salmon_flags"),
            (lambda r: r["design"].__setitem__("quantifier", "kallisto"), "quantifier"),
            (lambda r: r["salmon"].__setitem__("version", "1.10.3"), "Salmon '1.10.3'"),
            (lambda r: r.pop("genes_tpm"), "no genes_tpm rows"),
            (lambda r: r["design"].__setitem__("genes_simulated", ["CANX"]), "did not hold")):
        bad = copy.deepcopy(good)
        change(bad)
        with pytest.raises(gate.NotGradeable, match=re.escape(why)):
            gate.check_run(bad, PREDICTIONS)


def test_a_run_whose_numreads_rows_do_not_reproduce_the_stored_run_is_not_graded():
    run = _run()
    assert gate.check_reproduces(run) == 0.0
    run["genes"][0]["bias"] += gate.REPRODUCE_TOLERANCE / 2
    assert gate.check_reproduces(run) == pytest.approx(gate.REPRODUCE_TOLERANCE / 2)
    run["genes"][0]["bias"] += gate.REPRODUCE_TOLERANCE
    with pytest.raises(gate.NotGradeable, match="not made as stored"):
        gate.check_reproduces(run)
    run = _run()
    run["genes"][3]["n"] -= 1
    with pytest.raises(gate.NotGradeable, match="counts 29 replicates"):
        gate.check_reproduces(run)
    run = _run()
    run["genes"].pop()
    with pytest.raises(gate.NotGradeable, match="NumReads rows are for 48 genes"):
        gate.check_reproduces(run)


def test_the_command_exits_by_the_verdict(tmp_path, capsys):
    path = tmp_path / "run.json"
    path.write_text(json.dumps(_run()))
    assert gate.main([str(path)]) == 0
    out = capsys.readouterr().out
    assert "reproduce the stored uniform run (largest difference 0 log2)" in out
    assert "GATE PASSES: 3.0 makes `abundance` the default." in out
    path.write_text(json.dumps(_run(overrides={"SMN2": (0.5, 1.1)})))
    assert gate.main([str(path)]) == 1
    out = capsys.readouterr().out
    assert "The stop rule (step 4) stops the change: SMN2." in out
    assert "the control is not graded" in out and "GATE FAILS" in out
    run = _run()
    run["salmon"]["version"] = "2.6.0"
    path.write_text(json.dumps(run))
    assert gate.main([str(path)]) == 2
    assert "NOT GRADED" in capsys.readouterr().err
    assert hashlib.sha256((VALIDATION / "gate_3_0_predictions.json").read_bytes()).hexdigest() \
        in _text()
