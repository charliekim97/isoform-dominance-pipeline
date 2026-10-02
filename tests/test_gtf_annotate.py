"""``annotate --gtf``: groups proposed from a local GTF, with no network (issue #13).

``test_gtf_parity.py`` holds the GTF path to REST on a GENCODE 50 extract.  These are the
rules that extract cannot reach, the formats around it, and the command's failures.
"""
import gzip
import itertools
import json
import pathlib
import re

import pytest

from isoform_dominance import annotate, cli
from isoform_dominance import annotation_files as af

DATA = pathlib.Path(__file__).parent / "data" / "gencode_mini"
GTF = DATA / "gencode.v50.mini.gtf.gz"
REST = json.loads(gzip.decompress((DATA / "rest116_mini.json.gz").read_bytes()))
HEADER = "##description: test, version 50 (Ensembl 116)\n"


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Every test here runs with the Ensembl client unusable."""
    def refuse(*a, **k):
        raise AssertionError("a request was made")
    for name in ("get_json", "request_json", "resolve_server"):
        monkeypatch.setattr(annotate.ensembl, name, refuse)


def _cfg(gene, gtf=GTF, **kw):
    cfg = annotate.build_config_from_gtf(gene, str(gtf), **kw)
    return cfg, cfg.pop("annotation_source")


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def _gene(gid, name, chrom="chr1", strand="+", start=100, end=900, version=1):
    return ('%s\tHAVANA\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s.%d"; gene_type "protein_coding"; '
            'gene_name "%s";\n' % (chrom, start, end, strand, gid, version, name))


def _tx(gid, name, tid, exons, cds, chrom="chr1", strand="+", tags=(), version=1):
    """A transcript's lines; ``cds`` is ``[(start, end, frame)]``."""
    attrs = ('gene_id "%s.1"; transcript_id "%s.%d"; gene_type "protein_coding"; '
             'gene_name "%s"; transcript_type "protein_coding";%s'
             % (gid, tid, version, name, "".join(' tag "%s";' % t for t in tags)))
    lo, hi = min(s for s, _ in exons), max(e for _, e in exons)
    out = ["%s\tHAVANA\ttranscript\t%d\t%d\t.\t%s\t.\t%s\n" % (chrom, lo, hi, strand, attrs)]
    out += ["%s\tHAVANA\texon\t%d\t%d\t.\t%s\t.\t%s\n" % (chrom, s, e, strand, attrs)
            for s, e in exons]
    out += ["%s\tHAVANA\tCDS\t%d\t%d\t.\t%s\t%s\t%s\n" % (chrom, s, e, strand, f, attrs)
            for s, e, f in cds]
    return "".join(out)


# --------------------------------------------------------------------------- #
# the three rules
# --------------------------------------------------------------------------- #
def test_a_cds_that_starts_mid_codon_is_padded_and_a_partial_codon_dropped():
    """cds_start_NF: Ensembl pads the 5' end with (3 - frame) % 3 N's, then drops what is
    left of a codon at the 3' end."""
    assert af.protein_length([(1, 300, "0")], 1) == 100
    # 302 nt from frame 1: 2 N's, 304 nt, 101 codons and one base left over
    assert af.protein_length([(1, 102, "1"), (201, 400, "0")], 1) == 101
    assert af.protein_length([], 1) is None


def test_on_the_minus_strand_the_5_prime_cds_is_the_one_that_ends_furthest_along():
    # 203 nt; the 5'-most CDS on the minus strand is 500-600, in frame 2: one N
    cds = [(100, 201, "0"), (500, 600, "2")]
    assert af.protein_length(cds, -1) == 68
    assert af.protein_length(cds, 1) == 67           # read as plus strand: frame 0


def test_the_canonical_transcript_is_the_exact_tag_not_extended_and_not_mane(tmp_path):
    """GENCODE 50 added Ensembl_canonical_extended (14,432 transcripts), and 878 of its
    protein-coding genes have no MANE transcript."""
    g = "ENSG00000000003"
    gtf = _write(tmp_path, "c.gtf", HEADER + _gene(g, "CANON") + "".join(
        _tx(g, "CANON", "ENST0000000000%d" % i, [(100, 200), (300 + 100 * i, 900)],
            [(150, 200, "0"), (300 + 100 * i, 340 + 100 * i, "0")], tags=tags)
        for i, tags in enumerate([("Ensembl_canonical_extended",), ("MANE_Select",),
                                  ("basic", "Ensembl_canonical"), ("CCDS",)], 1)))
    rec, = af.scan(gtf, symbol="CANON")
    assert rec["canonical_transcript"] == "ENST00000000003.1"
    assert [t["is_canonical"] for t in rec["Transcript"]] == [0, 0, 1, 0]


def test_clusters_are_listed_by_content_whatever_order_the_transcripts_come_in():
    """REST at 116 and the GENCODE 50 GTF list a gene's transcripts in different orders
    (104 of 109 survey genes), and through 2.5 the order of clusters with as many
    transcripts followed it.  Both modes now list them by content."""
    txs = [{"id": "T1", "protein_aa": 900, "terminal_acceptor": 100, "is_canonical": True},
           {"id": "T2", "protein_aa": 700, "terminal_acceptor": 400, "is_canonical": False},
           {"id": "T3", "protein_aa": 800, "terminal_acceptor": 300, "is_canonical": False},
           {"id": "T4", "protein_aa": 800, "terminal_acceptor": 200, "is_canonical": False},
           {"id": "T5", "protein_aa": 600, "terminal_acceptor": 500, "is_canonical": False},
           {"id": "T6", "protein_aa": 600, "terminal_acceptor": 500, "is_canonical": False}]
    seen = set()
    for order in itertools.permutations(txs):
        clusters = annotate.cluster_by_terminal_exon({"transcripts": list(order)})
        seen.add(tuple((c["n"], c["rep_aa"], c["acceptor"]) for c in clusters))
    assert seen == {((2, 600, 500), (1, 900, 100), (1, 800, 200), (1, 800, 300),
                     (1, 700, 400))}


def test_tied_with_is_listed_by_content_too():
    txs = [{"id": "T%d" % i, "protein_aa": aa, "terminal_acceptor": acc,
            "is_canonical": i == 1}
           for i, (aa, acc) in enumerate([(900, 100), (800, 400), (800, 300), (800, 200)], 1)]
    for order in itertools.permutations(txs):
        cfg = annotate._config("X", "homo_sapiens", 116,
                               {"gene_id": "G", "transcripts": list(order)}, None)
        assert [c["terminal_acceptor"] for c in cfg["_proposal"]["tied_with"]] == [300, 400]
        assert [c["terminal_acceptor"] for c in cfg["_clusters"]] == [100, 200, 300, 400]


# --------------------------------------------------------------------------- #
# ids, regions and formats
# --------------------------------------------------------------------------- #
#: MT-ND4 as GENCODE writes it: one exon, and a CDS of 1,378 nt whose stop codon is
#: completed by polyadenylation; REST gives Translation.length 459.
MT_ND4 = (
    'chrM\tENSEMBL\tgene\t10760\t12137\t.\t+\t.\tgene_id "ENSG00000198886.2"; gene_type '
    '"protein_coding"; gene_name "MT-ND4"; level 3; hgnc_id "HGNC:7459";\n'
    + "".join(
        'chrM\tENSEMBL\t%s\t10760\t%d\t.\t+\t%s\tgene_id "ENSG00000198886.2"; transcript_id '
        '"ENST00000361381.2"; gene_type "protein_coding"; gene_name "MT-ND4"; transcript_type '
        '"protein_coding"; transcript_name "MT-ND4-201"; level 3; protein_id '
        '"ENSP00000354961.2"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select";\n'
        % (feature, end, frame)
        for feature, end, frame in (("transcript", 12137, "."), ("exon", 12137, "."),
                                    ("CDS", 12137, "0"), ("start_codon", 10762, "0"))))


def test_chrm_is_mt_and_mt_nd4_is_459_aa(tmp_path):
    gtf = _write(tmp_path, "mt.gtf", HEADER + MT_ND4)
    rec, = af.scan(gtf, symbol="MT-ND4")
    assert (rec["id"], rec["seq_region_name"]) == ("ENSG00000198886", "MT")
    assert rec["Transcript"][0]["Translation"] == {"length": 459}
    cfg, _ = _cfg("MT-ND4", gtf)
    assert cfg["groups"] == {"iso_459aa": ["ENST00000361381"]}
    assert "_gene_choice" not in cfg


def _ensembl_format(tmp_path, header="#!genome-build GRCh38.p14\n"
                                     "#!genebuild-last-updated 2026-03\n"):
    """The mini GTF as Ensembl writes its GTF: no ``chr``, versions in ``*_version``
    attributes, ``*_biotype`` for ``*_type``, and an ``#!`` header naming no release."""
    def unversion(m):
        return '%s "%s"; %s_version "%s";' % (m[1], m[2], m[1][:-3], m[3])
    out = [header]
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            p[0] = "MT" if p[0] == "chrM" else p[0][3:]
            p[8] = re.sub(r'\b(gene_id|transcript_id|exon_id|protein_id) "([^".]+)\.(\d+)";',
                          unversion, p[8])
            p[8] = p[8].replace("gene_type", "gene_biotype").replace("transcript_type",
                                                                     "transcript_biotype")
            out.append("\t".join(p) + "\n")
    return _write(tmp_path, "ensembl.gtf", "".join(out))


@pytest.mark.parametrize("gene", ["LEPR", "FOXO1", "CD99"])
def test_an_ensembl_gtf_proposes_what_the_gencode_one_does(tmp_path, gene):
    ens = _ensembl_format(tmp_path)
    notes = []
    got, src = _cfg(gene, ens, release=116, notes=notes)
    want, _ = _cfg(gene)
    assert got == want
    assert (src["gencode_release"], src["ensembl_release"], src["date"]) == (None, None,
                                                                            "2026-03")
    assert notes == ["%s names no Ensembl release in its header; the release is recorded as "
                     "116, from --ensembl-release" % ens]


def test_an_ensembl_gtf_gives_the_same_versioned_records(tmp_path):
    ens = _ensembl_format(tmp_path)
    for gid in ("ENSG00000116678", "ENSG00000002586"):
        a, = af.scan(str(GTF), gene_id=gid)
        b, = af.scan(ens, gene_id=gid)
        assert a == b


def test_a_header_with_no_release_records_none_and_says_so(tmp_path):
    ens = _ensembl_format(tmp_path, header="")
    notes = []
    cfg, src = _cfg("LEPR", ens, notes=notes)
    assert cfg["ensembl_release"] is None and src["ensembl_release"] is None
    assert len(notes) == 1 and "no release is recorded" in notes[0]


def test_the_header_gives_both_releases_and_the_provider():
    h = af.header(str(GTF))
    assert h == {"provider": "GENCODE", "date": "2026-04-08", "gencode_release": 50,
                 "ensembl_release": 116,
                 "description": "evidence-based annotation of the human genome (GRCh38), "
                                "version 50 (Ensembl 116)"}
    _, src = _cfg("LEPR")
    # every transcript of the gene, of every biotype, as REST's lookup gives them
    lepr = len(REST["lookups"]["ENSG00000116678"]["Transcript"])
    assert src == dict(kind="gtf", **af.provenance(str(GTF)), **h, n_transcripts=lepr)
    assert list(src) == ["kind", "file", "bytes", "sha256", "provider", "gencode_release",
                         "ensembl_release", "date", "description", "n_transcripts"]


def test_a_mouse_release_is_kept_as_written(tmp_path):
    gtf = _write(tmp_path, "m.gtf", "##description: evidence-based annotation of the mouse "
                 "genome (GRCm39), version M37 (Ensembl 114)\n")
    assert (af.header(gtf)["gencode_release"], af.header(gtf)["ensembl_release"]) == ("M37",
                                                                                     114)


# --------------------------------------------------------------------------- #
# reading the file
# --------------------------------------------------------------------------- #
def test_a_needle_that_straddles_a_block_boundary_is_found():
    """The GTF is read in blocks, and the part of a block after its last newline is
    carried into the next: a line, and a needle in it, is never cut in two."""
    raw = gzip.decompress(GTF.read_bytes())
    at = raw.index(b'gene_name "LEPR";')
    want = af.scan(str(GTF), symbol="LEPR")
    for block in (at + 5, at + 11, 97, 1000, 4099, 1 << 16):
        assert af.scan(str(GTF), symbol="LEPR", block=block) == want, block
    # the same through annotate
    assert _cfg("LEPR", block=at + 5) == _cfg("LEPR")


def test_the_last_line_is_read_without_a_final_newline(tmp_path):
    text = HEADER + MT_ND4
    gtf = _write(tmp_path, "nonl.gtf", text.rstrip("\n"))
    assert af.scan(gtf, symbol="MT-ND4", block=64) == af.scan(
        _write(tmp_path, "nl.gtf", text), symbol="MT-ND4")


def test_gzip_is_told_by_its_bytes_not_its_name(tmp_path):
    plain = tmp_path / "plain.gtf.gz"                 # not gzipped, whatever its name
    plain.write_bytes(gzip.decompress(GTF.read_bytes()))
    zipped = tmp_path / "zipped.txt"
    zipped.write_bytes(GTF.read_bytes())
    want, _ = _cfg("LEPR")
    assert _cfg("LEPR", plain)[0] == want and _cfg("LEPR", zipped)[0] == want


def test_a_truncated_gzip_is_one_line(tmp_path, capsys):
    cut = tmp_path / "cut.gtf.gz"
    cut.write_bytes(GTF.read_bytes()[:20000])
    rc = cli.main(["annotate", "--gene", "LEPR", "--gtf", str(cut), "--out",
                   str(tmp_path / "o.json")])
    err = capsys.readouterr().err
    assert rc == 1 and len(err.strip().splitlines()) == 1
    assert "is not a complete gzip file" in err


def test_a_symbol_is_found_ignoring_case_when_no_gene_has_it_exactly(tmp_path):
    """Ensembl's lookup/symbol ignores case; the exact name is looked for first, and the
    second, case-ignoring read of the file is said."""
    info = {}
    assert af.scan(str(GTF), symbol="LEPR", info=info) and info["case_insensitive"] is False
    assert af.scan(str(GTF), symbol="lepr", info=info) and info["case_insensitive"] is True
    notes = []
    got, _ = _cfg("lepr", notes=notes)
    want, _ = _cfg("LEPR")
    assert got == dict(want, gene="lepr")
    assert notes == ["%s has no gene named exactly lepr; LEPR was found ignoring case, as "
                     "Ensembl's REST lookup finds it" % GTF]


# --------------------------------------------------------------------------- #
# GENCODE basic
# --------------------------------------------------------------------------- #
def _basic(tmp_path):
    """The mini GTF as GENCODE's basic GTF holds it: the gene lines, and every line of a
    transcript tagged basic."""
    with gzip.open(GTF, "rt") as fh:
        keep = [x for x in fh if x.startswith("#") or "\tgene\t" in x
                or 'tag "basic";' in x]
    p = tmp_path / "gencode.v50.basic.annotation.gtf.gz"
    p.write_bytes(gzip.compress("".join(keep).encode()))
    return str(p)


def test_a_basic_gtf_is_refused(tmp_path, capsys):
    """Its header is the comprehensive one's; it leaves out transcripts the index holds,
    and at release 110 changed 51 of 109 proposals."""
    basic = _basic(tmp_path)
    with pytest.raises(af.BasicGTF, match="is a GENCODE basic GTF"):
        annotate.build_config_from_gtf("LEPR", basic)
    rc = cli.main(["annotate", "--gene", "LEPR", "--gtf", basic, "--out",
                   str(tmp_path / "o.json")])
    err = capsys.readouterr().err
    assert rc == 1 and len(err.strip().splitlines()) == 1
    assert "comprehensive gencode.vN.annotation.gtf.gz" in err
    assert not (tmp_path / "o.json").exists()


def _many(n, untagged=()):
    g = "ENSG00000000003"
    return HEADER + _gene(g, "ONE") + "".join(
        _tx(g, "ONE", "ENST%011d" % i, [(100, 900 - i)], [(100, 399, "0")],
            tags=() if i in untagged else ("basic",)) for i in range(1, n + 1))


def test_one_line_of_a_transcript_not_tagged_basic_makes_it_comprehensive(tmp_path):
    n = af.MIN_BASIC_TRANSCRIPTS
    with pytest.raises(af.BasicGTF, match="every one of its %d " % (3 * n)):
        af.scan(_write(tmp_path, "b.gtf", _many(n)), symbol="ONE")
    assert len(af.scan(_write(tmp_path, "c.gtf", _many(n, untagged={n})), symbol="ONE")) == 1


def test_an_extract_of_a_few_genes_is_not_called_basic(tmp_path):
    """Every single-transcript gene's transcript is basic, so a few genes cut from the
    comprehensive GTF can be all basic; GENCODE's basic GTF holds tens of thousands."""
    n = af.MIN_BASIC_TRANSCRIPTS - 1
    assert len(af.scan(_write(tmp_path, "few.gtf", _many(n)), symbol="ONE")) == 1


# --------------------------------------------------------------------------- #
# choosing the gene: the rule REST follows
# --------------------------------------------------------------------------- #
def test_a_symbol_whose_genes_are_all_off_the_reference_chromosomes_is_refused(tmp_path):
    g = "ENSG00000000004"
    text = HEADER + _gene(g, "ALTONLY", chrom="KI270713.1") + _tx(
        g, "ALTONLY", "ENST00000000004", [(100, 900)], [(100, 399, "0")], chrom="KI270713.1")
    gtf = _write(tmp_path, "alt.gtf", text)
    with pytest.raises(annotate.NotOnReference, match="ENSG00000000004 on KI270713.1"):
        annotate.build_config_from_gtf("ALTONLY", gtf)
    # --gene-id takes it anyway, and says where it is
    cfg, _ = _cfg("ALTONLY", gtf, gene_id=g)
    assert "is on KI270713.1, not a reference chromosome" in cfg["_gene_choice"]["reason"]
    # and the region rule is about GRCh38's names: not applied to another species
    cfg, _ = _cfg("ALTONLY", gtf, species="mus_musculus")
    assert cfg["gene_id"] == g


def test_gene_id_takes_one_gene_of_two(tmp_path):
    cfg, src = _cfg("HERC3", gene_id="ENSG00000287542.2")
    assert cfg["gene_id"] == "ENSG00000287542"
    assert cfg["_gene_choice"]["reason"] == "given by --gene-id"
    assert src["n_transcripts"] == cfg["_gene_choice"]["candidates"][0]["n_transcripts"]


def test_a_gene_id_the_gtf_lacks_or_of_another_gene_is_refused():
    with pytest.raises(ValueError, match="has no gene ENSG00000000009"):
        annotate.build_config_from_gtf("LEPR", str(GTF), gene_id="ENSG00000000009.1")
    with pytest.raises(ValueError, match="ENSG00000150907 is FOXO1, not LEPR"):
        annotate.build_config_from_gtf("LEPR", str(GTF), gene_id="ENSG00000150907")


def test_a_symbol_the_gtf_lacks_is_refused():
    with pytest.raises(ValueError, match="has no gene named NOSUCHGENE"):
        annotate.build_config_from_gtf("NOSUCHGENE", str(GTF))


def test_the_release_given_must_be_the_files():
    assert _cfg("LEPR", release=116)[0]["ensembl_release"] == 116
    with pytest.raises(ValueError, match=r"is Ensembl release 116 .*--ensembl-release 115 "
                                         r"contradicts it"):
        annotate.build_config_from_gtf("LEPR", str(GTF), release=115)


# --------------------------------------------------------------------------- #
# the command
# --------------------------------------------------------------------------- #
def _annotate(tmp_path, capsys, *argv):
    out = tmp_path / "cfg.json"
    rc = cli.main(["annotate", "--out", str(out), *argv])
    o, err = capsys.readouterr()
    return rc, (json.loads(out.read_text()) if out.exists() else None), err


def test_the_command_writes_the_config_with_no_request(tmp_path, capsys):
    rc, cfg, err = _annotate(tmp_path, capsys, "--gene", "CD99", "--gtf", str(GTF))
    assert rc == 0
    want, src = _cfg("CD99")
    assert cfg == dict(want, annotation_source=src)
    # the chrY gene and how to take it, as from REST
    assert "use --gene-id ENSG00000292348" in err


@pytest.mark.parametrize("argv, says", [
    (["--gene", "HERC3"], "HERC3 names 2 genes on the reference chromosomes: "
                          "ENSG00000138641 at 4:"),
    (["--gene", "LEPR", "--gene-id", "ENSG00000000009"], "has no gene ENSG00000000009"),
    (["--gene", "LEPR", "--gene-id", "ENSG00000150907"], "is FOXO1, not LEPR"),
    (["--gene", "LEPR", "--ensembl-release", "115"], "--ensembl-release 115 contradicts it"),
    (["--gene", "NOSUCHGENE"], "has no gene named NOSUCHGENE"),
])
def test_each_failure_is_one_line_and_exit_1(tmp_path, capsys, argv, says):
    rc, cfg, err = _annotate(tmp_path, capsys, "--gtf", str(GTF), *argv)
    assert rc == 1 and cfg is None
    assert len(err.strip().splitlines()) == 1 and says in err, err


def test_a_missing_gtf_is_one_line(tmp_path, capsys):
    rc, cfg, err = _annotate(tmp_path, capsys, "--gene", "LEPR", "--gtf",
                             str(tmp_path / "none.gtf.gz"))
    assert rc == 1 and len(err.strip().splitlines()) == 1 and "none.gtf.gz" in err


def test_the_notes_reach_stderr(tmp_path, capsys):
    rc, cfg, err = _annotate(tmp_path, capsys, "--gene", "lepr", "--gtf", str(GTF), "--json")
    assert rc == 0 and "NOTE: %s has no gene named exactly lepr" % GTF in err
    rc, cfg, err = _annotate(tmp_path, capsys, "--gene", "LEPR", "--gtf",
                             _ensembl_format(tmp_path, header=""))
    assert rc == 0 and cfg["ensembl_release"] is None
    assert "names no Ensembl release in its header; no release is recorded" in err
