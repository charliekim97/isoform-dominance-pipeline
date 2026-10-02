"""`annotate`: genes off the reference chromosomes, other species, odd symbols, and --json.

The reference-chromosome rule (:data:`annotate.GENE_RULE`) is a statement about GRCh38: a
human symbol whose genes all lie on alternate loci -- HLA-DRB3, GSTT1, KIR*, LILRA3 and
19 more protein-coding symbols in Ensembl 116's cdna.all -- was proposed from one of them,
silently, and the recommended reference-chromosome index does not contain it.  For another
species the human chromosome names were applied as they stand, so a zebrafish gene on
chr23 was dropped for a same-name gene on chr5.
"""
import http.client
import json
import urllib.request

import pytest

import test_gene_choice as gc
from isoform_dominance import annotate, cli, ensembl
from test_gene_choice import CD99_X, CD99_Y, Rest, _gene

# HLA-DRB3 at release 116: lookup/symbol gives a gene on an MHC alternate haplotype, and
# no gene of that name lies on a reference chromosome
DRB3_APD = _gene("ENSG00000231679", "HLA-DRB3", "CHR_HSCHR6_MHC_APD_CTG1", 3_946_000,
                 ["ENST00000000601", "ENST00000000602"])
DRB3_COX = _gene("ENSG00000230463", "HLA-DRB3", "CHR_HSCHR6_MHC_COX_CTG1", 3_985_000,
                 ["ENST00000000611", "ENST00000000612"])
# zebrafish: one name on chr23 and on chr5
DANIO_23 = _gene("ENSDARG00000000001", "lepr", "23", 1_000_000,
                 ["ENSDART00000000001", "ENSDART00000000002"])
DANIO_5 = _gene("ENSDARG00000000005", "lepr", "5", 2_000_000,
                ["ENSDART00000000005", "ENSDART00000000006"])


@pytest.fixture
def rest(monkeypatch):
    for g in (DRB3_APD, DRB3_COX, DANIO_23, DANIO_5):
        monkeypatch.setitem(gc.GENES, g["id"], g)

    def _serve(**kw):
        fake = Rest(**kw)
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        return fake
    return _serve


DRB3 = {"symbol_gives": {"HLA-DRB3": DRB3_APD["id"]},
        "xrefs": {"HLA-DRB3": [DRB3_APD["id"], DRB3_COX["id"]]}}


# --------------------------------------------------------------------------- #
# C1: a human symbol with no gene on a reference chromosome
# --------------------------------------------------------------------------- #
def test_a_symbol_only_on_alternate_loci_stops_with_the_candidates(rest, tmp_path, capsys):
    rest(**DRB3)
    with pytest.raises(annotate.NotOnReference) as e:
        annotate.build_config("HLA-DRB3")
    for g in (DRB3_APD, DRB3_COX):
        assert g["id"] in str(e.value) and g["seq_region_name"] in str(e.value)
    out = tmp_path / "c.json"
    assert cli.main(["annotate", "--gene", "HLA-DRB3", "--out", str(out)]) == 1
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and DRB3_APD["id"] in err[0] and "--gene-id" in err[0]
    assert not out.exists()


def test_even_when_xrefs_lists_only_the_one_lookup_gave(rest):
    rest(symbol_gives={"HLA-DRB3": DRB3_APD["id"]}, xrefs={"HLA-DRB3": [DRB3_APD["id"]]})
    with pytest.raises(annotate.NotOnReference):
        annotate.build_config("HLA-DRB3")


def test_gene_id_takes_one_of_them_and_says_the_index_lacks_it(rest, tmp_path, capsys):
    rest(**DRB3)
    cfg = annotate.build_config("HLA-DRB3", gene_id=DRB3_COX["id"])
    assert cfg["gene_id"] == DRB3_COX["id"]
    assert cli.main(["annotate", "--gene", "HLA-DRB3", "--gene-id", DRB3_COX["id"],
                     "--out", str(tmp_path / "c.json")]) == 0
    err = capsys.readouterr().err
    assert ("the recommended reference-chromosome index does not contain this gene" in err)


def test_a_reference_chromosome_gene_is_unchanged(rest):
    rest()
    cfg = annotate.build_config("CD99")
    assert cfg["gene_id"] == CD99_X["id"]                # the pair rule, as before
    assert annotate.build_config("LEPR")["gene_id"] == gc.LEPR["id"]


# --------------------------------------------------------------------------- #
# C2: the reference-chromosome rule is human; elsewhere two genes of a name are ambiguous
# --------------------------------------------------------------------------- #
def test_another_species_with_two_genes_of_a_name_is_ambiguous(rest):
    rest(symbol_gives={"lepr": DANIO_23["id"]},
         xrefs={"lepr": [DANIO_23["id"], DANIO_5["id"]]})
    with pytest.raises(annotate.AmbiguousGene) as e:
        annotate.build_config("lepr", species="danio_rerio")
    assert DANIO_23["id"] in str(e.value) and DANIO_5["id"] in str(e.value)
    assert "--gene-id" in str(e.value)


def test_another_species_with_one_gene_of_a_name_goes_on(rest):
    rest(symbol_gives={"lepr": DANIO_23["id"]}, xrefs={"lepr": [DANIO_23["id"]]})
    cfg = annotate.build_config("lepr", species="danio_rerio")
    assert cfg["gene_id"] == DANIO_23["id"]


# --------------------------------------------------------------------------- #
# C3: a symbol is quoted into the URL, and a URL that cannot be sent is not retried
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("symbol,quoted", [("HLA DRB3", "HLA%20DRB3"), ("A/B", "A%2FB")])
def test_a_symbol_is_quoted_into_the_url(monkeypatch, tmp_path, capsys, symbol, quoted):
    seen = []

    def fake(req, timeout=None, **_):
        seen.append(req.full_url)
        from test_ensembl_release import _http_error
        raise _http_error(req.full_url, 400, b'{"error":"No valid lookup found"}')
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    assert cli.main(["annotate", "--gene", symbol, "--out", str(tmp_path / "c.json")]) == 1
    assert seen == [ensembl.SERVER + "/lookup/symbol/homo_sapiens/%s?expand=1" % quoted]
    assert len(capsys.readouterr().err.strip().splitlines()) == 1


def test_an_invalid_url_is_not_retried(monkeypatch, tmp_path, capsys):
    calls = []

    def fake(req, timeout=None, **_):
        calls.append(req.full_url)
        raise http.client.InvalidURL("URL can't contain control characters")
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    monkeypatch.setattr("time.sleep", lambda s: None)
    assert cli.main(["annotate", "--gene", "LEPR", "--out", str(tmp_path / "c.json")]) == 1
    assert len(calls) == 1
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and "control characters" in err[0]
    assert "could not be sent" in err[0] and "network connection" not in err[0]


# --------------------------------------------------------------------------- #
# C4: --json keeps the notes, on stderr
# --------------------------------------------------------------------------- #
def test_json_still_gives_the_chrx_note_on_stderr(rest, tmp_path, capsys):
    rest()
    assert cli.main(["annotate", "--gene", "CD99", "--out", str(tmp_path / "c.json"),
                     "--json"]) == 0
    out = capsys.readouterr()
    assert json.loads(out.out)["gene_id"] == CD99_X["id"]
    assert "NOTE" in out.err and "--gene-id %s" % CD99_Y["id"] in out.err


def test_json_still_gives_the_tie_note_on_stderr(rest, monkeypatch, tmp_path, capsys):
    tied = _gene("ENSG00000000777", "TIED", "3", 10_000,
                 ["ENST00000000701", "ENST00000000702", "ENST00000000703"])
    tied["Transcript"][2]["Exon"][1] = {"start": 3000, "end": 3100}   # a third cluster
    monkeypatch.setitem(gc.GENES, tied["id"], tied)
    rest(symbol_gives={"TIED": tied["id"]}, xrefs={"TIED": [tied["id"]]})
    assert cli.main(["annotate", "--gene", "TIED", "--out", str(tmp_path / "c.json"),
                     "--json"]) == 0
    out = capsys.readouterr()
    assert json.loads(out.out)["_proposal"]["tied_with"]
    assert "chosen by a tie" in out.err
