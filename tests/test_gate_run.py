"""The 3.0 gate's run (``validation/data/sim/c1_d1_tpm.json``) and the two files computed to
describe it (``validation/data/README.md``, its last table).

CI runs no quantifier and has not the index's FASTA.  What it can hold: the stored run grades
as it graded when it was made; the per-transcript detail gives the run's own rows for its
genes, through the harness's own truth; and the whole-index file's gene-background half is
what the package gives now, through the call the harness made.  ``validation/whole_index.py``
is tried on a FASTA small enough to carry in a test.
"""
import gzip
import hashlib
import importlib.util
import json
import math
import pathlib
import statistics

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
DATA = VALIDATION / "data"


def _module(name):
    spec = importlib.util.spec_from_file_location("_" + name, VALIDATION / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


harness = _module("run_simulation")
gate = _module("gate_3_0")
summarize = _module("summarize")
whole_index = _module("whole_index")
RUN = json.loads((DATA / "sim" / "c1_d1_tpm.json").read_text(encoding="utf-8"))
DETAIL = json.loads((DATA / "sim" / "c1_d1_tpm_detail.json").read_text(encoding="utf-8"))
WHOLE = json.loads((DATA / "whole_index.json").read_text(encoding="utf-8"))
with gzip.open(DATA / "sweep116.json.gz", "rt", encoding="utf-8") as _f:
    SOURCE = json.load(_f)["genes"]


def _close(a, b):
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, rel_tol=1e-9)


# --------------------------------------------------------------------------- #
# the run
# --------------------------------------------------------------------------- #
def test_the_stored_run_grades_as_it_did(capsys):
    """Graded on 2026-10-23 with exit status 1: the stop rule, on SMN2."""
    assert gate.main([str(DATA / "sim" / "c1_d1_tpm.json")]) == 1
    out = capsys.readouterr().out
    assert "reproduce the stored uniform run (largest difference 0 log2)" in out
    assert "The stop rule (step 4) stops the change: SMN2." in out
    assert "The primary set fails; the control is not graded." in out
    assert out.rstrip().endswith("GATE FAILS: 3.0 keeps `fragment` as the default.")
    result = gate.grade(RUN, gate.frozen_predictions())
    assert result["stopped"] == ["SMN2"] and result["primary_sd_ok"] == 6
    smn2 = result["primary"]["SMN2"][0]
    assert (round(abs(smn2["bias"]), 3), round(smn2["realised_sd"], 3), smn2["n"]) == \
        (0.474, 0.751, 29)
    assert round(gate.frozen_predictions()["primary"]["SMN2"], 3) == 0.442
    assert [g for g, (_, b, s) in result["primary"].items() if not (b and s)] == ["SMN2"]
    assert result["control"] is None and result["passes"] is False


def test_the_runs_numreads_rows_measure_what_the_stored_uniform_runs_do():
    stored = {r["gene"]: r for r in json.loads(
        (DATA / "sim" / "c1_d1.json").read_text(encoding="utf-8"))["genes"]}
    rows = {r["gene"]: r for r in RUN["genes"]}
    assert set(rows) == set(stored) and len(rows) == 49
    for gene, r in rows.items():
        assert {k: r[k] for k in ("n", "realised_sd", "bias", "rmse")} == \
            {k: stored[gene][k] for k in ("n", "realised_sd", "bias", "rmse")}, gene


# --------------------------------------------------------------------------- #
# the per-transcript detail
# --------------------------------------------------------------------------- #
def _rows(gene):
    """The NumReads and TPM rows of ``gene`` recomputed from the detail, the TPM truth the
    harness's own (``run_simulation.tpm_truth``)."""
    g = DETAIL["genes"][gene]
    pc, groups = g["primary_comparison"], g["classes"]
    seqs = {t.split(".")[0]: s for t, s in SOURCE[gene]["sequences"].items()}
    per = {"genes": ([], []), "genes_tpm": ([], [])}
    for rep in g["replicates"]:
        q, f = rep["quant"], rep["fragments"]
        eff = {t: v[0] for t, v in q.items()}
        for key, col, truth in (("genes", 2, None), ("genes_tpm", 1, "tpm")):
            est = {c: sum(q[t][col] for t in groups[c] if t in q) for c in pc}
            if truth:
                tru, _ = harness.tpm_truth(groups, pc, seqs, f, eff, g["kept_copies"])
            else:
                tru = {c: sum(f[t] for t in groups[c]) for c in pc}
            if all(est[c] > 0 and tru[c] > 0 for c in pc):
                e = math.log2(est[pc[0]] / est[pc[1]])
                per[key][0].append(e)
                per[key][1].append(e - math.log2(tru[pc[0]] / tru[pc[1]]))
    return {key: {"n": len(e), "realised_sd": statistics.stdev(e),
                  "bias": statistics.median(x),
                  "rmse": math.sqrt(sum(v * v for v in x) / len(x))}
            for key, (e, x) in per.items()}


@pytest.mark.parametrize("gene", ["LDHA", "NDUFS1", "SMN2"])
def test_the_detail_gives_the_runs_rows(gene):
    assert DETAIL["run"] == {"file": "sim/c1_d1_tpm.json", "sha256": summarize.sha256_table(
        DATA / "README.md")["sim/c1_d1_tpm.json"]}
    assert len(DETAIL["genes"][gene]["replicates"]) == RUN["design"]["replicates"] == 30
    assert DETAIL["genes"][gene]["kept_copies"] == \
        RUN["effective_length_taken_from"].get(gene, {})
    mine = _rows(gene)
    for key in ("genes", "genes_tpm"):
        theirs = {r["gene"]: r for r in RUN[key]}[gene]
        assert mine[key]["n"] == theirs["n"], key
        for k in ("realised_sd", "bias", "rmse"):
            assert abs(mine[key][k] - theirs[k]) < 1e-12, (key, k)


def test_the_detail_has_every_replicates_reads_outside_the_49_genes():
    out = DETAIL["outside_num_reads"]
    assert len(out) == len(DETAIL["fragments_written"]) == 30 and all(r for r in out)
    assert all(n > 0 for r in out for n in r.values())
    # less than 1% of the fragments are given to records outside the 49 genes
    assert 0 < sum(sum(r.values()) for r in out) < 0.01 * sum(DETAIL["fragments_written"])
    for g in DETAIL["genes"].values():
        for rep in g["replicates"]:
            assert sum(rep["fragments"].values()) <= rep["gene_fragments_written"]
    # LDHA's and NDUFS1's copies outside the bundle go by those genes' symbols; SMN1, which
    # no fragment was drawn from, is not one of the 49
    assert {"LDHA", "NDUFS1", "SMN1"} <= set().union(*out)
    assert "SMN1" not in SOURCE


# --------------------------------------------------------------------------- #
# the whole-index file
# --------------------------------------------------------------------------- #
def test_the_whole_index_file_names_the_index_fasta_and_the_inputs():
    name, size, sha = whole_index.index_fasta()
    assert WHOLE["fasta"] == {"file": name, "bytes": size, "sha256": sha}
    assert sha == "683eb19310c40bf1396e4718f45afa2ce86755717c0990f47a171f535d248ea1"
    assert WHOLE["inputs"] == {"file": "sweep116.json.gz", "sha256": summarize.sha256_table(
        DATA / "README.md")["sweep116.json.gz"]}
    assert WHOLE["design"] == whole_index.DESIGN
    assert sorted(WHOLE["genes"]) == sorted(SOURCE)
    # the FASTA's records on chromosomes 1-22, X, Y and MT, 446,099 of its 465,769
    assert WHOLE["reference_chromosomes"] == {
        "records": 446099,
        "sha256": "50993c71f1ef9a5b27e2e492bcfba032721c3e784c545a74be6ba02aefa1720a"}


def test_against_the_whole_index_three_abundance_contrasts_are_not_estimable():
    """What the CHANGELOG says of the file: against the index's whole FASTA, SMN2's, CANX's
    and PTPRC's contrasts are not estimable as abundance and PTPRC's not as fragment
    totals; against its reference-chromosome records, only SMN2's, as abundance.  No
    contrast becomes estimable.  SMN1's transcripts are among SMN2's columns in both, and
    CANX's and PTPRC's own copies on scaffolds among theirs in the whole index."""
    genes = WHOLE["genes"]
    for b, est, lost, n_gene, n_after in (
            ("whole_index", "abundance", ["CANX", "PTPRC", "SMN2"], 48, 45),
            ("whole_index", "fragment", ["PTPRC"], 41, 40),
            ("reference_chromosomes", "abundance", ["SMN2"], 48, 47),
            ("reference_chromosomes", "fragment", [], 41, 41)):
        before = {g for g, v in genes.items() if v["gene_background"][est]["estimable"]}
        after = {g for g, v in genes.items() if v[b][est]["estimable"]}
        assert sorted(before - after) == lost and after <= before, (b, est)
        assert (len(before), len(after)) == (n_gene, n_after), (b, est)
    for b in ("whole_index", "reference_chromosomes"):
        smn2 = genes["SMN2"][b]["fasta_competitors"]
        assert sum(v["gene_symbol"] == "SMN1" for v in smn2.values()) >= 20, b
    for gene in ("CANX", "PTPRC"):
        comp = genes[gene]["whole_index"]["fasta_competitors"].values()
        own = [v for v in comp if v["gene_symbol"] == gene]
        assert own and all(v["location"].startswith("scaffold:") for v in own), gene
        assert not any(v["gene_symbol"] == gene for v in
                       genes[gene]["reference_chromosomes"]["fasta_competitors"].values())


def test_the_whole_index_files_gene_background_half_is_the_packages():
    """What ``whole_index.gene_background`` gives now, and the call is the harness's: the
    fragment and abundance figures equal ``run_simulation.predicted`` and
    ``predicted_abundance``."""
    for gene, entry in sorted(SOURCE.items()):
        now = whole_index.gene_background(entry)
        stored = WHOLE["genes"][gene]["gene_background"]
        for est in ("fragment", "abundance"):
            for k in ("verdict", "estimable", "beyond_linear"):
                assert stored[est][k] == now[est][k], (gene, est, k)
            for k in ("conditioning_factor", "min_resolvable_log2fc", "predicted_se"):
                v = now[est][k]
                assert _close(stored[est][k], v if math.isfinite(v) else None), (gene, est, k)
        frag = harness.predicted(SOURCE, gene, 10.0, 1.0, 1500.0)
        ab = harness.predicted_abundance(SOURCE, gene, 10.0, 1.0, 1500.0)
        assert now["fragment"]["predicted_se"] == frag["predicted_se"], gene
        assert now["abundance"]["predicted_se"] == ab["predicted_se"], gene
        assert now["fragment"]["verdict"] == frag["verdict"], gene


# a FASTA small enough for a test: half of one compared transcript on a scaffold, as an
# index of the whole assembly holds a gene's copies, half of the other under another gene's
# name, and records that share nothing -- on chromosomes 1, 22, X and MT, one without a
# symbol, and one on an unplaced scaffold
def _tiny_records(gene):
    seqs = {t.split(".")[0]: s for t, s in SOURCE[gene]["sequences"].items()}
    groups = SOURCE[gene]["config"]["groups"]
    pc = SOURCE[gene]["config"]["primary_comparison"]
    a, b = groups[pc[0]][0].split(".")[0], groups[pc[1]][0].split(".")[0]
    return [(">ENST90000000001.1 cdna scaffold:GRCh38:HSCHR_TEST:1:100:1 gene:ENSG9.1 "
             "gene_symbol:%s description:test" % gene, seqs[a][: len(seqs[a]) // 2]),
            (">ENST90000000002.1 cdna chromosome:GRCh38:7:1:100:1 gene:ENSG8.1 "
             "gene_symbol:PARALOG1 description:test", seqs[b][len(seqs[b]) // 2:]),
            (">ENST90000000003.1 cdna chromosome:GRCh38:1:1:100:1 gene:ENSG7.1 "
             "description:no symbol", "ACGT" * 40),
            (">ENST90000000004.1 cdna chromosome:GRCh38:22:1:100:1 gene:ENSG6.1 "
             "gene_symbol:FILLER22", "AAGC" * 40),
            (">ENST90000000005.1 cdna chromosome:GRCh38:X:1:100:1 gene:ENSG5.1 "
             "gene_symbol:FILLERX", "TTGA" * 40),
            (">ENST90000000006.1 cdna chromosome:GRCh38:MT:1:100:1 gene:ENSG4.1 "
             "gene_symbol:FILLERMT", "GGTC" * 40),
            (">ENST90000000007.1 cdna scaffold:GRCh38:KI270728.1:1:100:1 gene:ENSG3.1 "
             "gene_symbol:FILLERKI", "CATG" * 40)]


def _tiny(tmp_path, gene):
    path = tmp_path / "tiny.fa.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        for h, s in _tiny_records(gene):
            f.write("%s\n%s\n" % (h, s))
    return path


def _as_index(monkeypatch, path):
    data = path.read_bytes()
    monkeypatch.setattr(whole_index, "index_fasta",
                        lambda data_=None: ("tiny.fa.gz", len(data),
                                            hashlib.sha256(data).hexdigest()))


def test_whole_index_refuses_a_fasta_that_is_not_the_index(tmp_path):
    path = _tiny(tmp_path, "PTPRC")
    with pytest.raises(SystemExit, match="not the 683eb193"):
        whole_index.header(str(path))


def test_whole_index_gives_the_same_output_made_whole_or_in_parts(tmp_path, monkeypatch):
    path = _tiny(tmp_path, "CCND1")
    _as_index(monkeypatch, path)
    both = whole_index.run(str(path), ["CCND1", "COX5A"])
    first = tmp_path / "first.json"
    first.write_text(whole_index.as_json(whole_index.run(str(path), ["CCND1"])),
                     encoding="utf-8")
    parts = whole_index.run(str(path), ["CCND1", "COX5A"], [str(first)])
    assert whole_index.as_json(parts) == whole_index.as_json(both)
    comp = both["genes"]["CCND1"]["whole_index"]["fasta_competitors"]
    assert {r: (v["gene_symbol"], v["location"]) for r, v in comp.items()} == {
        "ENST90000000001": ("CCND1", "scaffold:GRCh38:HSCHR_TEST"),
        "ENST90000000002": ("PARALOG1", "chromosome:GRCh38:7")}
    assert all(v["windows"] > 0 for v in comp.values())
    # the reference-chromosome records leave out the two on scaffolds, and their SHA-256 is
    # that of the kept records' text
    kept = [(h, s) for h, s in _tiny_records("CCND1") if " chromosome:" in h]
    assert len(kept) == 5
    assert both["reference_chromosomes"] == {
        "records": 5,
        "sha256": hashlib.sha256("".join("%s\n%s\n" % r for r in kept).encode()).hexdigest()}
    assert set(both["genes"]["CCND1"]["reference_chromosomes"]["fasta_competitors"]) == {
        "ENST90000000002"}
    assert both["genes"]["COX5A"]["whole_index"]["fasta_competitors"] == {}
    other = tmp_path / "other.json"
    doc = json.loads(first.read_text(encoding="utf-8"))
    doc["isoform_dominance"] = "0.0"
    other.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(SystemExit, match="another FASTA, inputs, design or package"):
        whole_index.run(str(path), ["CCND1"], [str(other)])
