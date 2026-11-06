"""D1 (``validation/README.md``, "What produced SMN2's failure (D1)") as fixed before its run:
the records each arm leaves out, the checks of an arm's FASTA, index and ``quant.sf``, and the
grade.  CI has neither ``cdna.all`` nor Salmon, so the FASTA step is tried on a small FASTA and
the grade on runs made from the stored one.
"""
import copy
import gzip
import hashlib
import importlib.util
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"


def _module(name):
    spec = importlib.util.spec_from_file_location("_" + name, VALIDATION / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


d1 = _module("d1_smn2")
gate = _module("gate_3_0")
STORED = json.loads((VALIDATION / "data" / "sim" / "c1_d1_tpm.json").read_text(encoding="utf-8"))
README = (VALIDATION / "README.md").read_text(encoding="utf-8")
PREDICTIONS = gate.frozen_predictions()
SE = PREDICTIONS["primary"]["SMN2"]


def _section():
    return README[README.index("## What produced SMN2's failure (D1)"):]


# --------------------------------------------------------------------------- #
# what is fixed
# --------------------------------------------------------------------------- #
def test_the_left_out_lists_are_the_ones_fixed():
    for arm in d1.ARMS:
        ids = d1.LEFT_OUT[arm]
        assert len(set(ids)) == len(ids)
        assert d1._ids_sha256(ids) == d1.LEFT_OUT_SHA256[arm]
    assert set(d1.LEFT_OUT["ref"]) < set(d1.LEFT_OUT["gene"])
    assert (len(d1.LEFT_OUT["ref"]), len(d1.LEFT_OUT["gene"])) == (33, 56)


def test_the_readme_states_what_the_code_fixes():
    sec = _section()
    for arm in d1.ARMS:
        assert d1.LEFT_OUT_SHA256[arm] in sec, arm
        assert d1.ARM_FASTA_SHA256[arm] in sec, arm
        assert "{:,}".format(d1.NUM_REFS[arm]) in sec, arm
        assert "%d records" % len(d1.LEFT_OUT[arm]) in sec, arm
    assert "%.3f" % SE in sec
    assert d1.CDNA_SHA256 in (VALIDATION / "data" / "README.md").read_text(encoding="utf-8")


@pytest.mark.parametrize("header, ref, gene", [
    ("ENST1.1 cdna chromosome:GRCh38:5:1:2:1 gene:G1 gene_symbol:SMN1 description:x",
     False, True),
    ("ENST2.1 cdna scaffold:GRCh38:HSCHR5_1_CTG1_1:1:2:1 gene:G2 gene_symbol:SMN1",
     True, True),
    ("ENST3.1 cdna chromosome:GRCh38:5:1:2:1 gene:G3 gene_symbol:SMN2 description:x",
     False, False),
    ("ENST4.1 cdna scaffold:GRCh38:HSCHR5_2_CTG1_1:1:2:1 gene:G4 gene_symbol:SMN2",
     True, True),
    ("ENST5.1 cdna chromosome:GRCh38:9:1:2:1 gene:G5 gene_biotype:lncRNA", False, False),
    ("ENST6.1 cdna chromosome:GRCh38:5:1:2:1 gene:G6 gene_symbol:SMN1P1", False, False),
    ("ENST7.1 cdna chromosome:GRCh38:15:1:2:1 gene:G7 gene_symbol:X description:SMN1",
     False, False),
    ("ENST8.1 cdna chromosome:GRCh38:50:1:2:1 gene:G8 gene_symbol:SMN2", True, True),
])
def test_which_records_an_arm_leaves_out(header, ref, gene):
    assert d1.leaves_out(header, "ref") is ref
    assert d1.leaves_out(header, "gene") is gene


# --------------------------------------------------------------------------- #
# the arm's FASTA, index and quant.sf
# --------------------------------------------------------------------------- #
RECORDS = [
    ("ENST01.1 cdna chromosome:GRCh38:5:1:9:1 gene:A gene_symbol:SMN2", "ACGTACGTAA\nCCGG\n"),
    ("ENST02.2 cdna chromosome:GRCh38:5:1:9:1 gene:B gene_symbol:SMN1", "TTTTGGGGCC\n"),
    ("ENST03.1 cdna scaffold:GRCh38:HSCHR5_1_CTG1_1:1:9:1 gene:C gene_symbol:SMN2",
     "GGGGAAAACC\n"),
    ("ENST04.1 cdna scaffold:GRCh38:HSCHR5_1_CTG1_1:1:9:1 gene:D gene_symbol:SMN1",
     "CCCCAAAAGG\n"),
    ("ENST05.1 cdna chromosome:GRCh38:9:1:9:1 gene:E gene_symbol:OTHER", "ACACACACAC\n"),
]


def _write_cdna(path):
    with gzip.GzipFile(path, "wb", mtime=0) as f:
        for header, seq in RECORDS:
            f.write((">%s\n%s" % (header, seq)).encode("utf-8"))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fix(monkeypatch, cdna_sha):
    monkeypatch.setattr(d1, "CDNA_SHA256", cdna_sha)
    want = {}
    for arm in d1.ARMS:
        left = [h.split()[0].split(".")[0] for h, _ in RECORDS if d1.leaves_out(h, arm)]
        text = "".join(">%s\n%s" % (h, s) for h, s in RECORDS if not d1.leaves_out(h, arm))
        monkeypatch.setitem(d1.LEFT_OUT, arm, tuple(sorted(left)))
        monkeypatch.setitem(d1.LEFT_OUT_SHA256, arm, d1._ids_sha256(left))
        monkeypatch.setitem(d1.ARM_FASTA_SHA256, arm,
                            hashlib.sha256(text.encode("utf-8")).hexdigest())
        want[arm] = (sorted(left), text)
    return want


def test_an_arms_fasta_is_cdna_all_less_its_records(tmp_path, monkeypatch):
    want = _fix(monkeypatch, _write_cdna(tmp_path / "cdna.fa.gz"))
    assert want["ref"][0] == ["ENST03", "ENST04"]
    assert want["gene"][0] == ["ENST02", "ENST03", "ENST04"]
    for arm in d1.ARMS:
        out = tmp_path / ("%s.fa.gz" % arm)
        left = d1.write_fasta(str(tmp_path / "cdna.fa.gz"), arm, str(out))
        assert sorted(left) == want[arm][0]
        with gzip.open(out, "rt", encoding="utf-8", newline="") as f:
            assert f.read() == want[arm][1]
        assert not (tmp_path / ("%s.fa.gz.part" % arm)).exists()


@pytest.mark.parametrize("break_", ["input", "left out", "text"])
def test_the_fasta_step_refuses_another_input_or_result(tmp_path, monkeypatch, break_):
    _fix(monkeypatch, _write_cdna(tmp_path / "cdna.fa.gz"))
    if break_ == "input":
        monkeypatch.setattr(d1, "CDNA_SHA256", "0" * 64)
    elif break_ == "left out":
        monkeypatch.setitem(d1.LEFT_OUT_SHA256, "gene", d1._ids_sha256(["ENST03", "ENST04"]))
        monkeypatch.setitem(d1.LEFT_OUT, "gene", ("ENST03", "ENST04"))
    else:
        monkeypatch.setitem(d1.ARM_FASTA_SHA256, "gene", "0" * 64)
    out = tmp_path / "gene.fa.gz"
    with pytest.raises(d1.NotGradeable):
        d1.write_fasta(str(tmp_path / "cdna.fa.gz"), "gene", str(out))
    assert not out.exists() and not (tmp_path / "gene.fa.gz.part").exists()


def _info(tmp_path, **change):
    info = dict(d1.INDEX_INFO, num_refs=d1.NUM_REFS["ref"], seq_hash="x")
    info.update(change)
    (tmp_path / "info.json").write_text(json.dumps(info), encoding="utf-8")
    return str(tmp_path)


def test_an_arms_index_is_checked(tmp_path):
    assert d1.check_index("ref", _info(tmp_path))["num_refs"] == 453527
    for change in ({"num_refs": 453553}, {"num_refs": d1.NUM_REFS["gene"]},
                   {"keep_duplicates": True}, {"k": 25}, {"salmon_version": "1.10.3"},
                   {"num_decoys": 4}):
        with pytest.raises(d1.NotGradeable):
            d1.check_index("ref", _info(tmp_path, **change))


def _quant(tmp_path, names, gz):
    path = tmp_path / ("quant.sf.gz" if gz else "quant.sf")
    text = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n" + "".join(
        "%s\t100\t50\t1\t1\n" % n for n in names)
    if gz:
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return str(path)


@pytest.mark.parametrize("gz", [False, True])
def test_a_quant_sf_against_an_arms_index_is_checked(tmp_path, monkeypatch, gz):
    monkeypatch.setitem(d1.NUM_REFS, "gene", 3)
    assert d1.check_quant("gene", _quant(tmp_path, ["ENST9.1", "ENST8.2", "ENST7"], gz)) == 3
    with pytest.raises(d1.NotGradeable):
        d1.check_quant("gene", _quant(tmp_path, ["ENST9.1", "ENST8.2"], gz))
    left = d1.LEFT_OUT["gene"][0] + ".3"
    with pytest.raises(d1.NotGradeable):
        d1.check_quant("gene", _quant(tmp_path, ["ENST9.1", "ENST8.2", left], gz))


# --------------------------------------------------------------------------- #
# the grade
# --------------------------------------------------------------------------- #
def _run(bias=None, sd=None, drop=False):
    run = copy.deepcopy(STORED)
    rows = run["genes_tpm"]
    i = next(i for i, r in enumerate(rows) if r["gene"] == "SMN2")
    if drop:
        del rows[i]
    else:
        rows[i] = dict(rows[i], bias=bias, realised_sd=sd)
    return run


MEETS = (-0.20, 0.30)
MISSES = (-0.47, 0.75)


@pytest.mark.parametrize("ref, gene, exit_, words", [
    (MEETS, MEETS, 0, "met it."),
    (MISSES, MEETS, 0, "did not meet it."),
    (MEETS, MISSES, 1, "with ref it met it"),
    (MISSES, MISSES, 1, "with ref it did not meet it"),
])
def test_the_grade_is_gene_arms_and_reports_ref_beside_it(ref, gene, exit_, words):
    result = d1.grade({"ref": _run(*ref), "gene": _run(*gene)}, PREDICTIONS, STORED)
    assert result["exit"] == exit_
    text = d1.report(result)
    assert text.splitlines()[-1].startswith("D1 exit %d: %s" % (
        exit_, "THE INDEX." if exit_ == 0 else "NOT THE INDEX."))
    assert words in text.splitlines()[-1]
    # the stored run's own row, as graded on 2026-10-23, does not meet it
    assert "|bias| 0.474, SD 0.751, n 29" in text and not result["stored"]["meets"]
    assert result["others"]["ref"]["bias"][0] == 0.0


@pytest.mark.parametrize("bias, sd, met", [
    (-SE, 0.1, False), (SE, 0.1, False), (0.1, SE, False), (SE * 1.0000001, 0.1, False),
    (-SE * (1 - 1e-9), SE * (1 - 1e-9), True), (0.0, 0.0, True)])
def test_a_row_meets_the_predicted_se_only_below_it(bias, sd, met):
    assert d1.meets({"bias": bias, "realised_sd": sd, "n": 30}, SE) is met


def test_an_arm_without_smn2s_row_does_not_meet_it():
    result = d1.grade({"ref": _run(*MEETS), "gene": _run(drop=True)}, PREDICTIONS, STORED)
    assert result["exit"] == 1 and result["arms"]["gene"]["tpm"] is None
    assert "no row (fewer than 3 replicates count)" in d1.report(result)


@pytest.mark.parametrize("change", ["design", "salmon", "genes", "no tpm rows"])
def test_the_grade_refuses_runs_not_made_as_named(change):
    bad = _run(*MEETS)
    if change == "design":
        bad["design"]["bias_5p"] = 1.0
    elif change == "salmon":
        bad["salmon"]["version"] = "1.10.3"
    elif change == "genes":
        bad["design"]["genes_simulated"] = bad["design"]["genes_simulated"][:-1] + ["XYZ"]
    else:
        del bad["genes_tpm"]
    with pytest.raises(d1.NotGradeable):
        d1.grade({"ref": _run(*MEETS), "gene": bad}, PREDICTIONS, STORED)


def test_main_grades_and_names_its_exit_status(tmp_path, capsys):
    paths = {}
    for name, run in (("meets", _run(*MEETS)), ("misses", _run(*MISSES)),
                      ("bad", dict(_run(*MEETS), salmon={"version": "0"}))):
        paths[name] = tmp_path / ("%s.json" % name)
        paths[name].write_text(json.dumps(run), encoding="utf-8")
    assert d1.main(["grade", str(paths["misses"]), str(paths["meets"])]) == 0
    assert capsys.readouterr().out.rstrip().splitlines()[-1].startswith("D1 exit 0: THE INDEX.")
    assert d1.main(["grade", str(paths["meets"]), str(paths["misses"])]) == 1
    assert "D1 exit 1: NOT THE INDEX." in capsys.readouterr().out
    assert d1.main(["grade", str(paths["meets"]), str(paths["bad"])]) == 2
    assert capsys.readouterr().out.startswith("D1 exit 2: NOT GRADED.")


def test_a_refused_check_says_so_on_stderr_and_exits_2(tmp_path, capsys):
    info = dict(d1.INDEX_INFO, num_refs=1)
    (tmp_path / "info.json").write_text(json.dumps(info), encoding="utf-8")
    assert d1.main(["index", "ref", str(tmp_path)]) == 2
    out = capsys.readouterr()
    assert out.out == "" and out.err.startswith("d1_smn2.py index refused: ")
    info["num_refs"] = d1.NUM_REFS["ref"]
    (tmp_path / "info.json").write_text(json.dumps(info), encoding="utf-8")
    assert d1.main(["index", "ref", str(tmp_path)]) == 0
    assert capsys.readouterr().out.startswith("ref: index ")


# --------------------------------------------------------------------------- #
# D1's result
# --------------------------------------------------------------------------- #
#: D1's section of validation/README.md, from its title to its result, as fixed before its
#: run (ef2c4fe): the protocol is not edited after it
D1_TEXT_SHA256 = "25bdc9c56d057584611a3319880ce31b4447b3a2166dda0955c79f7d9f13f932"
SIM = VALIDATION / "data" / "sim"
RUN_PATHS = {arm: SIM / ("d1_%s_tpm.json" % arm) for arm in d1.ARMS}
RUNS = {arm: json.loads(path.read_text(encoding="utf-8")) for arm, path in RUN_PATHS.items()}
RESULT = d1.grade(RUNS, PREDICTIONS, STORED)


def _flat(text):
    return " ".join(text.split())


def _result():
    sec = README[README.index("### D1's result"):]
    return sec[:sec.index("\n## ")]


def test_the_protocol_is_the_text_fixed_before_its_run():
    t = README[README.index("## What produced SMN2's failure (D1)"):
               README.index("\n### D1's result")]
    assert t.isascii() and hashlib.sha256(t.encode("ascii")).hexdigest() == D1_TEXT_SHA256


def test_the_stored_runs_grade_as_the_readme_says(capsys):
    assert RESULT["exit"] == 0
    assert RESULT["arms"]["ref"]["meets"] and not RESULT["stored"]["meets"]
    assert d1.main(["grade", str(RUN_PATHS["ref"]), str(RUN_PATHS["gene"])]) == 0
    last = capsys.readouterr().out.rstrip().splitlines()[-1]
    assert last.startswith("D1 exit 0: THE INDEX.") and last.endswith("(ref) it met it.")
    assert "graded `data/sim/d1_ref_tpm.json` and `data/sim/d1_gene_tpm.json` with exit status 0" \
        in _flat(_result())


def _cells(entry):
    tpm, num = entry["tpm"], entry["num_reads"]
    return " | ".join(["%.3f" % abs(tpm["bias"]), "%.3f" % tpm["realised_sd"], "%d" % tpm["n"],
                       "yes" if entry["meets"] else "no",
                       "%.3f" % abs(num["bias"]), "%.3f" % num["realised_sd"]])


def test_the_readmes_table_is_the_grade():
    sec = _result()
    assert "| below %.3f |" % SE in sec
    assert "| the gate's run | none | %s |" % _cells(RESULT["stored"]) in sec
    for arm in d1.ARMS:
        line = "| `%s` | %d | %s |" % (arm, len(d1.LEFT_OUT[arm]), _cells(RESULT["arms"][arm]))
        assert line in sec, line


def test_the_readme_says_what_ref_kept_of_smn1():
    """The records `gene` leaves out and `ref` keeps are SMN1's on chromosome 5: the columns
    that make SMN2's contrast not estimable against the reference chromosomes are among them,
    and the others stand, under SMN1's names, for SMN2's transcripts of the same sequence."""
    kept = set(d1.LEFT_OUT["gene"]) - set(d1.LEFT_OUT["ref"])
    whole = json.loads((VALIDATION / "data" / "whole_index.json").read_text(encoding="utf-8"))
    cols = whole["genes"]["SMN2"]["reference_chromosomes"]["fasta_competitors"]
    smn1 = {x for x, c in cols.items() if c["gene_symbol"] == "SMN1"}
    assert smn1 <= kept
    assert all(cols[x]["location"] == "chromosome:GRCh38:5" for x in smn1)
    detail = json.loads((SIM / "c1_d1_tpm_detail.json").read_text(encoding="utf-8"))
    stand_ins = {x: s for entry in detail["outside_stand_ins"] for x, s in entry.items()}
    rest = kept - smn1
    assert rest and all(stand_ins[x]["symbol"] == "SMN1" for x in rest)
    assert all(v["gene"] == d1.GENE for x in rest for v in stand_ins[x]["stands_for"].values())
    sec = _flat(_result())
    assert ("SMN1's %d transcripts on chromosome 5 stayed in that index, %d of them among the "
            "columns that make SMN2's contrast not estimable" % (len(kept), len(smn1))) in sec
    assert ("so those %d were columns the quantifier could leave empty; the other %d stand, "
            "under SMN1's names, for %d of SMN2's transcripts of the same sequence"
            % (len(smn1), len(rest), len(rest))) in sec
    assert "SMN1" not in STORED["design"]["genes_simulated"]


def test_the_readme_draws_what_the_grade_says():
    """Exit status 0 is "the failure was the index's" (the protocol's words), and `ref`
    meeting the figure is what the README's sentences on it say."""
    assert RESULT["exit"] == 0 and RESULT["arms"]["ref"]["meets"]
    sec = _flat(_result())
    assert sec.count("**The failure was the index's.**") == 1
    assert ("With only the two genes' copies on alternate loci left out, records the README's "
            "\"Index scope\" also leaves out, it met it too") in sec
    shows = _flat(README[README.index("**What this shows, and what it does not.**"):
                         README.index("## What produced SMN2's failure (D1)")])
    assert ("D1, below, measured it: with only the index's copies of SMN1 and SMN2 on alternate "
            "loci left out, SMN1's transcripts on chromosome 5 still in it and no fragment drawn "
            "from SMN1, SMN2's summed TPM met its predicted SE.") in shows


def test_the_data_readme_names_each_arms_index():
    text = _flat((VALIDATION / "data" / "README.md").read_text(encoding="utf-8"))
    assert ("its index without the %d SMN1 and SMN2 records off chromosome 5 (%s records)"
            % (len(d1.LEFT_OUT["ref"]), "{:,}".format(d1.NUM_REFS["ref"]))) in text
    assert ("its index without every SMN1 record and SMN2's records off chromosome 5, %d (%s "
            "records)" % (len(d1.LEFT_OUT["gene"]), "{:,}".format(d1.NUM_REFS["gene"]))) in text
    for arm in d1.ARMS:
        row = [line for line in text.split("| `sim/") if line.startswith("d1_%s_tpm.json" % arm)]
        assert len(row) == 1, arm


def test_the_other_genes_figures_are_the_grades():
    sec = _flat(_result())
    o = RESULT["others"]
    assert ("in `ref` the largest bias change is %.4f (%s) and the largest SD change %.4f (%s); "
            "in `gene`, %.4f (%s) and %.4f (%s)"
            % (o["ref"]["bias"] + o["ref"]["realised_sd"]
               + o["gene"]["bias"] + o["gene"]["realised_sd"])) in sec
    n = {name: {r["gene"]: r["n"] for r in run["genes_tpm"]}
         for name, run in [("stored", STORED)] + sorted(RUNS.items())}
    moved = [g for g in n["stored"] if g != d1.GENE
             and not n["stored"][g] == n["ref"][g] == n["gene"][g]]
    assert moved == [o["gene"]["realised_sd"][1]]
    g = moved[0]
    assert n["stored"][g] == n["ref"][g]
    assert ("%s's row counts %d replicates in `gene` against %d in the gate's run and in `ref`"
            % (g, n["gene"][g], n["stored"][g])) in sec
    # every other row moved; the NumReads rows' largest bias change is the same in either arm
    worst = set()
    for arm in d1.ARMS:
        for key in ("genes", "genes_tpm"):
            want = {r["gene"]: r for r in STORED[key] if r["gene"] != d1.GENE}
            got = {r["gene"]: r for r in RUNS[arm][key] if r["gene"] != d1.GENE}
            assert all(got[x]["bias"] != want[x]["bias"] for x in want), (arm, key)
        want = {r["gene"]: r for r in STORED["genes"] if r["gene"] != d1.GENE}
        got = {r["gene"]: r for r in RUNS[arm]["genes"] if r["gene"] != d1.GENE}
        d, x = max((abs(got[x]["bias"] - want[x]["bias"]), x) for x in want)
        worst.add(("%.3f" % d, x))
    assert len(worst) == 1
    ((d, x),) = worst
    assert "The other 48 genes' rows, against the gate's run, all moved, by little." in sec
    assert "`NumReads`: the largest bias change is %s (%s) in either arm." % (d, x) in sec


def test_the_gates_numbers_in_either_arm_are_the_readmes():
    seen = set()
    for arm in d1.ARMS:
        g = gate.grade(RUNS[arm], PREDICTIONS)
        bias_ok = sum(1 for _, b, _ in g["primary"].values() if b)
        seen.add((bias_ok, len(g["primary"]), g["primary_sd_ok"], tuple(g["stopped"]),
                  g["control"]["bias_ok"], g["control"]["sd_ok"]))
    assert len(seen) == 1
    ((bias_ok, n_primary, sd_ok, stopped, c_bias, c_sd),) = seen
    assert stopped == ()
    assert ("the primary set's |bias| below its predicted SE for %d of %d and its SD for %d, and "
            "the control's for %d and %d, where the gate asks for all %d, %d, %d and %d"
            % (bias_ok, n_primary, sd_ok, c_bias, c_sd, n_primary, gate.PRIMARY_SD_AT_LEAST,
               gate.CONTROL_BIAS_AT_LEAST, gate.CONTROL_SD_AT_LEAST)) in _flat(_result())


def test_the_top_level_readme_quotes_the_grade():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    text = _flat(text.replace("\n> ", "\n"))
    s, r = RESULT["stored"]["tpm"], RESULT["arms"]["ref"]["tpm"]
    assert RESULT["arms"]["ref"]["meets"] and not RESULT["stored"]["meets"]
    assert "SMN1" not in STORED["design"]["genes_simulated"]
    assert ("In the simulation of `validation/`, which draws no fragment from SMN1, SMN2's class "
            "ratio of summed TPM missed its predicted SE with SMN1's and SMN2's copies on "
            "alternate loci in the index (|bias| %.3f log2, SD %.3f, against %.3f) and met it "
            "with them left out (%.3f and %.3f; D1)"
            % (abs(s["bias"]), s["realised_sd"], SE, abs(r["bias"]), r["realised_sd"])) in text
    assert ("with only SMN1's and SMN2's copies on alternate loci left out, records \"Index "
            "scope\" below also leaves out, SMN2's summed TPM met its predicted SE, in a "
            "simulation that draws no fragment from SMN1") in text
    assert ("in that simulation, which draws no fragment from SMN1, leaving out the index's "
            "copies of SMN1 and SMN2 on alternate loci was enough for that row to meet its "
            "predicted SE, SMN1's transcripts on chromosome 5 still in the index") in text


def test_the_paper_gives_d1s_result():
    """paper.md: SMN2's figures held without the index's copies of SMN1 and SMN2 on scaffolds,
    the records `ref` leaves out; SMN1 had no simulated reads."""
    paper = _flat((ROOT / "paper.md").read_text(encoding="utf-8"))
    assert ("A check fixed before its run saw SMN2's figures hold without the index's copies of "
            "SMN1 and SMN2 on scaffolds, which the package's documentation excludes; the "
            "paralogue SMN1 had no simulated reads.") in paper
    assert RESULT["exit"] == 0 and RESULT["arms"]["ref"]["meets"]
    assert not RESULT["stored"]["meets"]
    assert "SMN1" not in STORED["design"]["genes_simulated"]
    whole = json.loads((VALIDATION / "data" / "whole_index.json").read_text(encoding="utf-8"))
    cols = whole["genes"]["SMN2"]["whole_index"]["fasta_competitors"]
    copies = {x for x, c in cols.items() if c["gene_symbol"] in ("SMN1", "SMN2")
              and c["location"].startswith("scaffold:")}
    assert copies and copies <= set(d1.LEFT_OUT["ref"])


def test_the_cli_note_gives_d1s_result():
    """The gene-background NOTE when the abundance estimand decides: after SMN2's failure,
    what `ref` gave, and that no fragment was drawn from SMN1."""
    from isoform_dominance import cli
    note = cli._WHOLE_INDEX_MEASURED["abundance"]
    assert RESULT["arms"]["ref"]["meets"] and not RESULT["stored"]["meets"]
    assert "SMN1" not in STORED["design"]["genes_simulated"]
    assert ("SMN2's is the contrast whose summed TPM failed that gate; with SMN1's and SMN2's "
            "copies on alternate loci left out of the index, its summed TPM met its predicted "
            "SE, in a simulation that draws no fragment from SMN1 (README, ") in note
