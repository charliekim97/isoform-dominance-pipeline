"""The GTF path gives what REST gives, on a GENCODE 50 extract (issue #13).

Ensembl 116 (June 2026) is the last release the legacy platform publishes, and its REST
API is kept "for e116 for long term use" only; the new platform has none.  A GENCODE GTF
and its transcript FASTA are what is left for later releases, and for every release
offline.  The fixtures (``tests/data/gencode_mini/README.md``) are all of release 116 /
GENCODE 50:

  gencode.v50.mini.gtf.gz              every GTF line of LEPR, FOXO1, STK11, AXIN1, GSK3B,
                                       CD99 (chrX and chrY genes) and HERC3 (two genes)
  gencode.v50.mini.transcripts.fa.gz   their records from gencode.v50.transcripts.fa.gz
  rest116_mini.json.gz                 the REST lookups of the same genes, trimmed to the
                                       fields the package reads, and MD5 of REST cDNA
"""
import gzip
import hashlib
import json
import pathlib

import pytest

from isoform_dominance import annotate
from isoform_dominance import annotation_files as af

DATA = pathlib.Path(__file__).parent / "data" / "gencode_mini"
GTF = DATA / "gencode.v50.mini.gtf.gz"
FASTA = DATA / "gencode.v50.mini.transcripts.fa.gz"
REST = json.loads(gzip.decompress((DATA / "rest116_mini.json.gz").read_bytes()))
BY_NAME = {}
for _g in REST["lookups"].values():
    BY_NAME.setdefault(_g["display_name"], []).append(_g)


@pytest.fixture
def rest(monkeypatch):
    """annotate's REST calls answered from the recorded lookups."""
    def get(path, **kw):
        if path.startswith("/info/data"):
            return {"releases": [REST["release"]]}
        sym = path.split("?")[0].split("/")[4] if "/symbol/" in path else None
        if path.startswith("/lookup/symbol/"):
            gid = REST["lookup_symbol"].get(sym) or BY_NAME[sym][0]["id"]
            return REST["lookups"][gid]
        if path.startswith("/xrefs/symbol/"):
            return [{"type": "gene", "id": g["id"]} for g in BY_NAME.get(sym, [])]
        if path.startswith("/lookup/id/"):
            return REST["lookups"][path.split("/")[3].split("?")[0]]
        raise AssertionError(path)
    monkeypatch.setattr(annotate, "_get", get)
    monkeypatch.setattr(annotate, "_post", lambda path, body, **kw:
                        {i: REST["lookups"][i] for i in body["ids"] if i in REST["lookups"]})
    monkeypatch.setattr(annotate.ensembl, "resolve_server", lambda release=None, **kw: None)


def _file_cfg(gene, **kw):
    cfg = annotate.build_config_from_gtf(gene, str(GTF), **kw)
    src = cfg.pop("annotation_source")
    assert (src["gencode_release"], src["ensembl_release"]) == (50, 116)
    return cfg


@pytest.mark.parametrize("gene", ["LEPR", "FOXO1", "STK11", "AXIN1", "GSK3B", "CD99"])
def test_the_file_path_proposes_exactly_what_rest_proposes(rest, gene):
    assert _file_cfg(gene) == annotate.build_config(gene)


def test_of_a_par_pair_the_chrx_gene_is_proposed_from(rest):
    """Parity alone cannot see this: both paths share the X/Y rule."""
    cfg = _file_cfg("CD99")
    assert cfg["gene_id"] == "ENSG00000002586"
    assert cfg["_gene_choice"]["chosen"] == "ENSG00000002586"
    assert [c["gene_id"] for c in cfg["_gene_choice"]["candidates"]] == [
        "ENSG00000002586", "ENSG00000292348"]


@pytest.mark.parametrize("gid", sorted(REST["lookups"]))
def test_every_field_annotate_reads_is_the_same(gid):
    """Per transcript, not per proposal: a canonical flag or a protein length can be wrong
    without moving any proposal in this extract (Ensembl_canonical_extended, cds_start_NF)."""
    ref = REST["lookups"][gid]
    got = next(g for g in af.scan(str(GTF), gene_id=gid) if g["id"] == gid)
    for key in ("display_name", "seq_region_name", "start", "end", "strand",
                "canonical_transcript"):
        assert got[key] == ref[key], key
    a = annotate.transcripts_of(ref, ref["display_name"])
    b = annotate.transcripts_of(got, ref["display_name"])
    assert sorted(a["transcripts"], key=lambda t: t["id"]) \
        == sorted(b["transcripts"], key=lambda t: t["id"])
    assert {t["id"].split(".")[0] for t in ref["Transcript"]} \
        == {t["id"] for t in got["Transcript"]}


def test_the_cdna_is_rests_byte_for_byte():
    checked = 0
    for g in af.scan(str(GTF), symbol="LEPR") + af.scan(str(GTF), gene_id="ENSG00000002586"):
        for t, s in af.sequences_for(g, str(FASTA)).items():
            if t in REST["cdna_md5"]:
                assert hashlib.md5(s.encode()).hexdigest() == REST["cdna_md5"][t], t
                checked += 1
    assert checked > 50


def test_a_symbol_of_two_unrelated_genes_stops_with_both():
    with pytest.raises(annotate.AmbiguousGene, match="ENSG00000138641.*ENSG00000287542"):
        annotate.build_config_from_gtf("HERC3", str(GTF))


def _par_gtf(tmp_path, with_y):
    lines = []
    for chrom, suffix in (("chrX", ""), ("chrY", "_PAR_Y"))[:2 if with_y else 1]:
        g = ('gene_id "ENSG00000000001.5%s"; gene_type "protein_coding"; gene_name "PARG1";'
             % suffix)
        lines.append("%s\tHAVANA\tgene\t100\t900\t.\t+\t.\t%s\n" % (chrom, g))
        for tid, acc in (("ENST00000000001", 700), ("ENST00000000002", 800)):
            t = '%s transcript_id "%s.1%s"; transcript_type "protein_coding";%s' % (
                g, tid, suffix, ' tag "Ensembl_canonical";' if tid.endswith("1") else "")
            lines.append("%s\tHAVANA\ttranscript\t100\t900\t.\t+\t.\t%s\n" % (chrom, t))
            for s_, e in ((100, 200), (acc, 900)):
                lines.append("%s\tHAVANA\texon\t%d\t%d\t.\t+\t.\t%s\n" % (chrom, s_, e, t))
            lines.append("%s\tHAVANA\tCDS\t150\t200\t.\t+\t0\t%s\n" % (chrom, t))
            lines.append("%s\tHAVANA\tCDS\t%d\t%d\t.\t+\t0\t%s\n" % (chrom, acc, acc + 32, t))
    p = tmp_path / ("par_%s.gtf" % with_y)
    p.write_text("##description: test, version 43 (Ensembl 109)\n" + "".join(lines))
    cfg = annotate.build_config_from_gtf("PARG1", str(p))
    cfg.pop("annotation_source")
    return cfg


def test_pre_44_par_y_records_change_nothing(tmp_path):
    """GENCODE 25-43 give the chrY PAR copy the chrX ids plus _PAR_Y (GENCODE _README.TXT);
    REST before 110 has no chrY gene for it, so the config is the chrX-only one."""
    assert _par_gtf(tmp_path, True) == _par_gtf(tmp_path, False)
