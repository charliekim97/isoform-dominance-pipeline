"""Tests for five mutants that the 2.4.0 suite let live (post-release audit, appendix A).

Each test here fails under the mutant named beside it; see the 2.4.1 round's mutation
table.  M27 (``extract`` in donor order) is killed by ``test_extract_order.py``.
"""
import json
import urllib.request

import pytest

from isoform_dominance import annotate, cli, io
from test_gene_choice import CD99_X, CD99_Y, Rest, _gene
import test_gene_choice as gc
from test_identical_copies import SEQS, T1, TWIN, _case, _run


# M2: the configured side of the identical-sequence test compared as given, not upper-cased
def test_a_lowercase_configured_sequence_is_identical_to_its_uppercase_record(tmp_path,
                                                                             capsys):
    lower = {t: s.lower() for t, s in SEQS.items()}
    _, res, err = _run(_case(tmp_path, {TWIN: SEQS[T1].upper()}, seqs=lower), capsys)
    assert res["background"]["identical_to_configured"] == {TWIN: T1}
    assert "1 record(s) with the sequence of a configured transcript" in err


# M4: the pseudoautosomal rule taken for any set of genes whose regions are X and Y
def test_two_genes_on_x_and_one_on_y_are_ambiguous_not_a_par_pair(monkeypatch):
    x2 = _gene("ENSG00000002587", "CD99", "X", 3_000_000, ["ENST00000000105"])
    monkeypatch.setitem(gc.GENES, x2["id"], x2)
    monkeypatch.setattr(urllib.request, "urlopen", Rest(
        xrefs={"CD99": [CD99_X["id"], CD99_Y["id"], x2["id"]]}))
    with pytest.raises(annotate.AmbiguousGene) as e:
        annotate.build_config("CD99")
    for g in (CD99_X, CD99_Y, x2):
        assert g["id"] in str(e.value)


# M36: an offline save takes the config's gene_id when nothing was fetched
def test_an_offline_save_records_the_config_gene_id(tmp_path, capsys, monkeypatch):
    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    argv = _case(tmp_path)
    cfg = json.loads(open(argv[2]).read())
    cfg["gene_id"] = "ENSG00000000077"
    open(argv[2], "w").write(json.dumps(cfg))
    saved = tmp_path / "s.json"
    assert cli.main(argv + ["--save-inputs", str(saved)]) == 0
    assert json.loads(saved.read_text())["gene_id"] == "ENSG00000000077"


# M37: a rerun from saved inputs that saves again keeps the saved gene_id
def test_a_resave_keeps_the_gene_id_of_the_file_it_ran_from(tmp_path, capsys, monkeypatch):
    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    argv = _case(tmp_path)                       # its config records no gene_id
    first, again = tmp_path / "first.json", tmp_path / "again.json"
    io.save_inputs(str(first), {"sequences": dict(SEQS), "background_sequences": {},
                                "gene_id": "ENSG00000000077", "k": 31, "window": 31,
                                "canonical": True}, json.loads(open(argv[2]).read()), 116, "t")
    rerun = [a for a in argv if a not in ("--sequences", argv[argv.index("--sequences") + 1])]
    assert cli.main(rerun + ["--inputs", str(first), "--save-inputs", str(again)]) == 0
    assert json.loads(again.read_text())["gene_id"] == "ENSG00000000077"
