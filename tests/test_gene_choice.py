"""Which gene a symbol means, when the reference chromosomes hold more than one of that name.

Ensembl REST ``lookup/symbol`` returns one gene per symbol.  For the pseudoautosomal genes
CD99, CRLF2, CSF2RA, IL3RA and SHOX it returned the chrY copy (releases 110 and 116), and
for HERC3 and DUSP13B the newer of two genes of that name on one chromosome.  `annotate`
took whatever it was given.  A pseudoautosomal gene's chrY copy has the chrX copy's
sequence, and Salmon keeps only the first of identical sequences: in a GENCODE-built index
that is chrX, so a config proposed from the chrY gene named no transcript quant.sf has.

The HTTP layer is replaced at ``urllib.request.urlopen``, so these tests see every request
actually made.
"""
import json
import random
import urllib.request

import pytest

from isoform_dominance import annotate, cli

SERVER = "https://rest.ensembl.org"


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _gene(gid, name, region, start, tids):
    """A gene payload as ``lookup/...?expand=1`` gives it: two terminal-exon clusters."""
    tx = []
    for i, t in enumerate(tids):
        acc = 1000 if i % 2 == 0 else 2000
        tx.append({"id": t + ".1", "biotype": "protein_coding", "is_canonical": int(i == 0),
                   "Translation": {"length": 300 if acc == 1000 else 250},
                   "Exon": [{"start": 1, "end": 100}, {"start": acc, "end": acc + 100}]})
    return {"id": gid, "display_name": name, "seq_region_name": region, "start": start,
            "end": start + 50_000, "strand": 1, "canonical_transcript": tids[0] + ".1",
            "Transcript": tx}


CD99_X = _gene("ENSG00000002586", "CD99", "X", 2_690_988,
               ["ENST00000000101", "ENST00000000102", "ENST00000000103", "ENST00000000104"])
CD99_Y = _gene("ENSG00000292348", "CD99", "Y", 2_690_988,
               ["ENST00000000201", "ENST00000000202", "ENST00000000203", "ENST00000000204"])
HERC3_OLD = _gene("ENSG00000138641", "HERC3", "4", 88_592_372,
                  ["ENST00000000301", "ENST00000000302", "ENST00000000303", "ENST00000000304"])
HERC3_NEW = _gene("ENSG00000287542", "HERC3", "4", 88_523_810, ["ENST00000000311"])
LEPR = _gene("ENSG00000116678", "LEPR", "1", 65_420_543,
             ["ENST00000000401", "ENST00000000402", "ENST00000000403"])
SMN1 = _gene("ENSG00000172062", "SMN1", "5", 70_925_030,
             ["ENST00000000501", "ENST00000000502", "ENST00000000503"])
SMN2 = _gene("ENSG00000205571", "SMN2", "5", 70_049_635, ["ENST00000000511"])
SMN1_ALT = _gene("ENSG00000275349", "SMN1", "HSCHR5_1_CTG1_1", 473_491,
                 ["ENST00000000521"])
GENES = {g["id"]: g for g in (CD99_X, CD99_Y, HERC3_OLD, HERC3_NEW, LEPR, SMN1, SMN2,
                              SMN1_ALT)}
# what lookup/symbol gave on 2026-09-29, and what xrefs/symbol lists (SMN1's has SMN2's)
SYMBOL_GIVES = {"CD99": CD99_Y["id"], "HERC3": HERC3_NEW["id"], "LEPR": LEPR["id"],
                "SMN1": SMN1["id"]}
XREFS = {"CD99": [CD99_X["id"], CD99_Y["id"]], "HERC3": [HERC3_OLD["id"], HERC3_NEW["id"]],
         "LEPR": [LEPR["id"], "LRG_1234"], "SMN1": [SMN1["id"], SMN2["id"], SMN1_ALT["id"]]}

# a pseudoautosomal pair has one sequence: the chrY transcript is the chrX one's
SEQS = {}
for n, t in enumerate(t for g in GENES.values() for t in g["Transcript"]):
    SEQS[t["id"].split(".")[0]] = _seq(900 + 150 * (n % 2), 100 + n)
for x, y in zip(CD99_X["Transcript"], CD99_Y["Transcript"], strict=True):
    SEQS[y["id"].split(".")[0]] = SEQS[x["id"].split(".")[0]]


class Rest:
    """lookup/symbol, xrefs/symbol, lookup/id (GET and POST), cDNA and the release; every
    request is logged as (method, path without the query, body)."""

    def __init__(self, symbol_gives=SYMBOL_GIVES, xrefs=XREFS):
        self.symbol_gives, self.xrefs = symbol_gives, xrefs
        self.calls = []

    def __call__(self, req, timeout=None, **_):
        from test_ensembl_http import _Resp
        path = req.full_url[len(SERVER):]
        body = json.loads(req.data) if req.data else None
        self.calls.append((req.get_method(), path.split("?")[0], body))
        parts = path.split("?")[0].split("/")
        if path.startswith("/info/data"):
            return _Resp(json.dumps({"releases": [116]}))
        if path.startswith("/lookup/symbol/"):
            return _Resp(json.dumps(GENES[self.symbol_gives[parts[4]]]))
        if path.startswith("/xrefs/symbol/"):
            return _Resp(json.dumps([{"type": "gene", "id": i}
                                     for i in self.xrefs.get(parts[4], [])]))
        if path.startswith("/lookup/id/"):
            return _Resp(json.dumps(GENES[parts[3].split(".")[0]]))
        if path.startswith("/lookup/id") and req.get_method() == "POST":
            return _Resp(json.dumps({i: GENES[i] for i in body["ids"] if i in GENES}))
        if path.startswith("/sequence/id") and req.get_method() == "POST":
            return _Resp(json.dumps([{"query": i, "id": i, "seq": SEQS[i.split(".")[0]]}
                                     for i in body["ids"]]))
        raise AssertionError("unexpected request %s %s" % (req.get_method(), path))

    def paths(self):
        return [(m, p) for m, p, _ in self.calls]


@pytest.fixture
def rest(monkeypatch):
    def _serve(**kw):
        fake = Rest(**kw)
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        return fake
    return _serve


def _ids(gene):
    return {t["id"].split(".")[0] for t in gene["Transcript"]}


def _configured(cfg):
    return {t for ids in cfg["groups"].values() for t in ids}


# --------------------------------------------------------------------------- #
# annotate
# --------------------------------------------------------------------------- #
def test_a_pseudoautosomal_pair_uses_the_chrx_gene(rest):
    rest()
    cfg = annotate.build_config("CD99")
    assert cfg["gene_id"] == CD99_X["id"]
    assert _configured(cfg) <= _ids(CD99_X)
    choice = cfg["_gene_choice"]
    assert choice["chosen"] == CD99_X["id"]
    assert [c["gene_id"] for c in choice["candidates"]] == [CD99_X["id"], CD99_Y["id"]]
    assert "pseudoautosomal" in choice["reason"] and "GENCODE" in choice["reason"]


def test_the_cli_names_the_chry_gene_and_how_to_take_it(rest, tmp_path, capsys):
    rest()
    assert cli.main(["annotate", "--gene", "CD99", "--out", str(tmp_path / "c.json")]) == 0
    err = capsys.readouterr().err
    assert CD99_X["id"] in err and CD99_Y["id"] in err
    assert "--gene-id %s" % CD99_Y["id"] in err


def test_same_name_genes_that_are_not_a_pair_stop_with_the_candidates(rest, tmp_path, capsys):
    rest()
    with pytest.raises(annotate.AmbiguousGene) as e:
        annotate.build_config("HERC3")
    msg = str(e.value)
    for g in (HERC3_OLD, HERC3_NEW):
        assert g["id"] in msg and "4:%d" % g["start"] in msg
    assert "4 transcripts" in msg and "1 transcript" in msg
    assert "--gene-id" in msg
    assert cli.main(["annotate", "--gene", "HERC3", "--out", str(tmp_path / "h.json")]) == 1
    err = capsys.readouterr().err
    assert HERC3_OLD["id"] in err and HERC3_NEW["id"] in err


def test_gene_id_chooses_and_asks_nothing_else(rest, tmp_path):
    fake = rest()
    cfg = annotate.build_config("HERC3", gene_id=HERC3_OLD["id"])
    assert cfg["gene_id"] == HERC3_OLD["id"]
    assert _configured(cfg) <= _ids(HERC3_OLD)
    assert cfg["_gene_choice"]["reason"] == "given by --gene-id"
    assert [p for m, p in fake.paths() if not p.startswith("/info")] \
        == ["/lookup/id/%s" % HERC3_OLD["id"]]
    assert cli.main(["annotate", "--gene", "HERC3", "--gene-id", HERC3_OLD["id"],
                     "--out", str(tmp_path / "h.json")]) == 0


def test_a_gene_id_of_another_gene_is_refused(rest):
    rest()
    with pytest.raises(ValueError, match="is CD99, not LEPR"):
        annotate.build_config("LEPR", gene_id=CD99_X["id"])


def test_one_candidate_costs_one_request_more_and_no_lookup_id(rest):
    """xrefs/symbol is the only request that lists every gene of a name, so it is made on
    every run: one request more than 2.3 made (lookup/symbol, /info/data).  lookup/id is
    asked only when xrefs/symbol lists a gene besides the one lookup/symbol gave."""
    fake = rest()
    cfg = annotate.build_config("LEPR")
    assert fake.paths() == [("GET", "/lookup/symbol/homo_sapiens/LEPR"),
                            ("GET", "/xrefs/symbol/homo_sapiens/LEPR"),
                            ("GET", "/info/data")]
    assert cfg["gene_id"] == LEPR["id"]
    assert "_gene_choice" not in cfg


def test_another_name_or_a_non_reference_region_is_not_a_candidate(rest):
    """xrefs/symbol/SMN1 lists SMN2 and SMN1's alternate-locus copy as well."""
    fake = rest()
    cfg = annotate.build_config("SMN1")
    assert cfg["gene_id"] == SMN1["id"]
    assert "_gene_choice" not in cfg
    assert ("POST", "/lookup/id") in fake.paths()


def test_an_empty_xrefs_answer_goes_on_with_what_lookup_gave_and_says_so(
        rest, tmp_path, capsys):
    rest(xrefs={})
    cfg = annotate.build_config("CD99")
    assert cfg["gene_id"] == CD99_Y["id"]
    assert "not looked for" in cfg["_gene_choice"]["reason"]
    assert cli.main(["annotate", "--gene", "CD99", "--out", str(tmp_path / "c.json")]) == 0
    assert "not looked for" in capsys.readouterr().err
