"""Index scope: same-name copies on non-reference regions, and which index quantified whom.

GENCODE's transcripts.fa gained scaffold, patch and alternate-locus transcripts at
release 48, and Ensembl's cdna.all has them too (release 116 checked).  A gene with a copy
there is in the FASTA twice, under its own name and another gene id; an index built from
it lets the quantifier split reads between the two, and `extract` counts only the
configured ids.  A same-name gene on a reference chromosome is not such a copy, but only
an Ensembl header says where a gene lies.
"""
import json
import os
import random

import pytest

from isoform_dominance import cli, extract, identifiability, index_scope


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


_SHARED = _seq(1500, 31)
_ALT_A = _seq(800, 32)
_ALT_B = _seq(1000, 33)
_OTHER = _seq(900, 34)
SEQS = {"ENST00000000001": _SHARED + _ALT_A, "ENST00000000002": _SHARED + _ALT_B}
GROUPS = {"A": ["ENST00000000001"], "B": ["ENST00000000002"]}


def _gencode(tid, gid, name, seq):
    return ">%s.1|%s.3|-|-|%s-201|%s|%d|protein_coding|\n%s\n" % (
        tid, gid, name, name, len(seq), seq)


def _ensembl(tid, gid, name, region, seq, kind="chromosome"):
    return (">%s.1 cdna %s:GRCh38:%s:100:200:1 gene:%s.3 gene_biotype:protein_coding "
            "transcript_biotype:protein_coding gene_symbol:%s description:x\n%s\n"
            % (tid, kind, region, gid, name, seq))


# the reference gene, ENSG00000000010, and an unrelated gene
_REF = (_gencode("ENST00000000001", "ENSG00000000010", "GENEX", SEQS["ENST00000000001"])
        + _gencode("ENST00000000002", "ENSG00000000010", "GENEX", SEQS["ENST00000000002"])
        + _gencode("ENST00000000009", "ENSG00000000090", "OTHER", _OTHER))
# ... and a copy of it on an alternate locus: same name, another gene id
_COPY = _gencode("ENST00000000005", "ENSG00000000050", "GENEX", _SHARED + _ALT_A[:400])
# ... and, instead, the chrY copy of a pseudoautosomal gene: same name, same sequence, its
# own gene id, as GENCODE 44, 48 and 50 write CD99 or SHOX -- on a reference chromosome
_PAR_Y = _gencode("ENST00000000006", "ENSG00000000060", "GENEX", SEQS["ENST00000000001"])


def _case(tmp_path, fasta_text, gene="GENEX"):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"gene": gene, "groups": GROUPS,
                               "primary_comparison": ["A", "B"]}))
    sq = tmp_path / "seqs.json"
    sq.write_text(json.dumps(SEQS))
    fa = tmp_path / "bg.fa"
    fa.write_text(fasta_text)
    return ["identifiability", "--config", str(cfg), "--sequences", str(sq),
            "--no-gene-background", "--background-fasta", str(fa)]


# ---- header parsing -------------------------------------------------------- #
def test_a_gencode_header_gives_transcript_gene_id_and_name():
    rec = index_scope.parse_header(
        "ENST00000380707.9|ENSG00000172062.19|OTTHUMG00000099361.6|"
        "OTTHUMT00000216881.2|SMN1-201|SMN1|1570|protein_coding|")
    assert rec == {"transcript": "ENST00000380707", "gene_id": "ENSG00000172062",
                   "gene_name": "SMN1", "region": None}


def test_an_ensembl_header_gives_the_region_too():
    rec = index_scope.parse_header(
        "ENST00000631131.2 cdna scaffold:GRCh38:HSCHR5_1_CTG1_1:473491:500558:1 "
        "gene:ENSG00000275349.4 gene_biotype:protein_coding "
        "transcript_biotype:protein_coding gene_symbol:SMN1 description:survival of "
        "motor neuron 1, telomeric [Source:HGNC Symbol;Acc:HGNC:11117]")
    assert rec == {"transcript": "ENST00000631131", "gene_id": "ENSG00000275349",
                   "gene_name": "SMN1", "region": "HSCHR5_1_CTG1_1"}


@pytest.mark.parametrize("head", ["DECOY", "chr1", "ENST00000000001.1", "a|b|c", "",
                                  "ENST1 cdna chromosome:GRCh38:1:1:2:1"])
def test_a_header_in_neither_format_is_not_parsed(head):
    assert index_scope.parse_header(head) is None


def test_a_par_y_transcript_is_the_same_gene_not_a_copy():
    """Older GENCODE releases (not 44, 48 or 50) wrote the chrY copy of a pseudoautosomal
    gene as <id>.<v>_PAR_Y; that is the gene itself, and must not be reported as a copy."""
    recs = [index_scope.parse_header(h) for h in (
        "ENST00000381192.10|ENSG00000182378.15|-|-|PLCXD1-201|PLCXD1|3163|protein_coding|",
        "ENST00000381187.8|ENSG00000182378.15|-|-|PLCXD1-202|PLCXD1|1801|protein_coding|",
        # an unconfigured transcript's PAR_Y record: only the gene id says whose it is
        "ENST00000381187.8_PAR_Y|ENSG00000182378.15_PAR_Y|-|-|PLCXD1-202|PLCXD1|1801|"
        "protein_coding|")]
    assert index_scope.same_name_copies(recs, ["ENST00000381192"]) == []


# ---- identifiability --background-fasta ------------------------------------ #
def test_a_gencode_fasta_with_a_same_name_gene_warns_that_it_may_be_a_copy(
        tmp_path, capsys):
    rc = cli.main(_case(tmp_path, _REF + _COPY) + ["--json"])
    out, err = capsys.readouterr()
    copies = json.loads(out)["background"]["same_name_copies"]
    assert copies == [{"gene_id": "ENSG00000000050", "gene_name": "GENEX", "region": None,
                       "transcripts": ["ENST00000000005"]}]
    assert rc != 1          # a warning, not an error
    assert "WARNING" in err
    assert "ENSG00000000050" in err
    assert "may be copies" in err
    assert "reference-chromosome" in err


def test_a_gencode_par_pair_warns_with_the_caveat(tmp_path, capsys):
    """A GENCODE header does not say where a transcript lies, so the chrY copy of a
    pseudoautosomal gene cannot be told from a copy on an alternate locus. It is reported,
    as a possible copy, and the warning says what else it may be."""
    cli.main(_case(tmp_path, _REF + _PAR_Y) + ["--json"])
    out, err = capsys.readouterr()
    assert json.loads(out)["background"]["same_name_copies"] == [
        {"gene_id": "ENSG00000000060", "gene_name": "GENEX", "region": None,
         "transcripts": ["ENST00000000006"]}]
    assert "WARNING" in err
    assert "may be copies" in err
    assert "pseudoautosomal" in err
    assert "not an index-scope problem" in err


def _ensembl_fasta(other_region, kind="chromosome"):
    """GENEX on chromosome X, and another gene id named GENEX on ``other_region``."""
    return (_ensembl("ENST00000000001", "ENSG00000000010", "GENEX", "X",
                     SEQS["ENST00000000001"])
            + _ensembl("ENST00000000002", "ENSG00000000010", "GENEX", "X",
                       SEQS["ENST00000000002"])
            + _ensembl("ENST00000000005", "ENSG00000000050", "GENEX", other_region,
                       _SHARED, kind=kind))


@pytest.mark.parametrize("region", ["HSCHR1_1_CTG3", "HG2290_PATCH", "KI270711.1"])
def test_an_ensembl_fasta_with_a_copy_names_its_region(region, tmp_path, capsys):
    cli.main(_case(tmp_path, _ensembl_fasta(region, kind="scaffold")) + ["--json"])
    out, err = capsys.readouterr()
    assert json.loads(out)["background"]["same_name_copies"] == [
        {"gene_id": "ENSG00000000050", "gene_name": "GENEX", "region": region,
         "transcripts": ["ENST00000000005"]}]
    assert "WARNING" in err and region in err
    # the region is known, so there is nothing to hedge
    assert "may be copies" not in err
    assert "pseudoautosomal" not in err


@pytest.mark.parametrize("region", ["Y", "X", "10", "MT"])
def test_an_ensembl_same_name_gene_on_a_reference_chromosome_is_silent(
        region, tmp_path, capsys):
    """The chrY copy of a pseudoautosomal gene (CD99, SHOX) and a distinct gene sharing the
    name (HERC3 on 4, DUSP13B on 10) are on reference chromosomes: an index built from
    reference chromosomes keeps them, so they are not what the check is for."""
    cli.main(_case(tmp_path, _ensembl_fasta(region)) + ["--json"])
    out, err = capsys.readouterr()
    assert json.loads(out)["background"]["same_name_copies"] == []
    assert "WARNING" not in err


def test_a_clean_fasta_is_silent(tmp_path, capsys):
    cli.main(_case(tmp_path, _REF) + ["--json"])
    out, err = capsys.readouterr()
    assert json.loads(out)["background"]["same_name_copies"] == []
    assert "WARNING" not in err


def test_a_fasta_in_neither_format_is_silent_and_runs(tmp_path, capsys):
    fa = ">DECOY\n%s\n>chr1 some text\n%s\n>x|y\n%s\n" % (_OTHER, _OTHER, _OTHER)
    rc = cli.main(_case(tmp_path, fa) + ["--json"])
    out, err = capsys.readouterr()
    assert rc == cli.EXIT_OK
    assert json.loads(out)["background"]["same_name_copies"] == []
    assert "WARNING" not in err


def test_the_copy_check_is_by_the_configured_transcripts_gene_not_by_the_config_name(
        tmp_path, capsys):
    """A hand-written config may name the gene differently, or not at all; the name that
    matters is the one the FASTA gives the configured transcripts."""
    cli.main(_case(tmp_path, _REF + _COPY, gene=None) + ["--json"])
    out, _ = capsys.readouterr()
    assert [c["gene_id"] for c in json.loads(out)["background"]["same_name_copies"]] \
        == ["ENSG00000000050"]


def test_no_background_fasta_reports_no_copies():
    res = identifiability.analyze({"groups": GROUPS, "primary_comparison": ["A", "B"]},
                                  sequences=SEQS, background_gene_transcripts=False)
    assert res["background"]["same_name_copies"] == []


# ---- extract: which index quantified each donor ---------------------------- #
def _quantdir(tmp_path, donors, names=None):
    """``donors``: {donor: meta_info dict, or None for no aux_info/meta_info.json}."""
    qd = tmp_path / "quant"
    names = names or ["ENST00000000001.1", "ENST00000000002.1"]
    for d, meta in donors.items():
        (qd / d).mkdir(parents=True)
        (qd / d / "quant.sf").write_text(
            "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
            + "".join("%s\t1000\t800\t%d\t10\n" % (n, 5 + i) for i, n in enumerate(names)))
        if meta is not None:
            (qd / d / "aux_info").mkdir()
            (qd / d / "aux_info" / "meta_info.json").write_text(json.dumps(meta))
    sm = tmp_path / "sm.csv"
    sm.write_text("donor,condition\n" + "".join("%s,control\n" % d for d in donors))
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"gene": "GENEX", "groups": GROUPS,
                               "primary_comparison": ["A", "B"]}))
    return ["extract", "--config", str(cfg), "--quantdir", str(qd), "--samplemap", str(sm),
            "--cohort", "C", "--out", str(tmp_path / "pd.csv")]


def _meta(seq_hash, targets=642692):
    return {"salmon_version": "1.10.3", "index_seq_hash": seq_hash,
            "index_name_hash": "n" + seq_hash, "num_valid_targets": targets,
            "keep_duplicates": False, "num_processed": 1}


def test_extract_records_the_index_in_a_sidecar(tmp_path):
    argv = _quantdir(tmp_path, {"D1": _meta("58f7"), "D2": _meta("58f7")})
    assert cli.main(argv) == cli.EXIT_OK
    side = json.load(open(str(tmp_path / "pd.csv") + ".index.json"))
    assert side["format"] == index_scope.INDEX_FORMAT
    assert side["index_seq_hashes"] == ["58f7"]
    assert side["mixed"] is False
    assert side["missing_meta_info"] == []
    assert side["donors"]["D1"] == {"salmon_version": "1.10.3", "index_seq_hash": "58f7",
                                    "index_name_hash": "n58f7",
                                    "num_valid_targets": 642692, "keep_duplicates": False,
                                    "index_decoy_seq_hash": None, "num_decoy_targets": None}
    assert side["unreadable_meta_info"] == {}
    # the per-donor CSV itself is untouched: no comment lines for `stats` to trip on
    assert open(tmp_path / "pd.csv").read().splitlines()[0].startswith("cohort,donor,")


def test_extract_stops_when_donors_were_quantified_against_different_indexes(
        tmp_path, capsys):
    argv = _quantdir(tmp_path, {"D1": _meta("58f7"), "D2": _meta("fdc2", 654828)})
    assert cli.main(argv) == 1
    err = capsys.readouterr().err
    assert "58f7" in err and "fdc2" in err
    assert "D1" in err and "D2" in err
    assert "--allow-mixed-index" in err
    assert not os.path.exists(tmp_path / "pd.csv")


def test_allow_mixed_index_passes_and_says_so_in_the_sidecar(tmp_path, capsys):
    argv = _quantdir(tmp_path, {"D1": _meta("58f7"), "D2": _meta("fdc2", 654828)})
    assert cli.main(argv + ["--allow-mixed-index"]) == cli.EXIT_OK
    assert "WARNING" in capsys.readouterr().err
    side = json.load(open(str(tmp_path / "pd.csv") + ".index.json"))
    assert side["mixed"] is True
    assert side["index_seq_hashes"] == ["58f7", "fdc2"]


def test_the_library_call_refuses_a_mixed_cohort_too(tmp_path):
    _quantdir(tmp_path, {"D1": _meta("58f7"), "D2": _meta("fdc2")})
    cfg = json.load(open(tmp_path / "cfg.json"))
    with pytest.raises(index_scope.MixedIndexError):
        extract.run(cfg, str(tmp_path / "quant"), str(tmp_path / "sm.csv"), "C",
                    str(tmp_path / "pd.csv"))


def test_no_meta_info_is_a_warning_not_an_error(tmp_path, capsys):
    argv = _quantdir(tmp_path, {"D1": None, "D2": _meta("58f7")})
    assert cli.main(argv) == cli.EXIT_OK
    err = capsys.readouterr().err
    assert "WARNING" in err and "D1" in err and "meta_info.json" in err
    side = json.load(open(str(tmp_path / "pd.csv") + ".index.json"))
    assert side["missing_meta_info"] == ["D1"]
    assert side["donors"]["D1"] is None
    assert side["mixed"] is False


def test_extract_warns_on_a_copy_named_in_full_gencode_headers(tmp_path, capsys):
    """An index built without --gencode keeps the whole header as the target name."""
    names = [line[1:] for line in (_REF + _COPY).splitlines() if line.startswith(">")]
    argv = _quantdir(tmp_path, {"D1": _meta("fdc2")}, names=names)
    assert cli.main(argv) == cli.EXIT_OK
    err = capsys.readouterr().err
    assert "WARNING" in err and "ENSG00000000050" in err
    assert "may be copies" in err and "pseudoautosomal" in err
    # the configured transcripts are still found under their full-header names
    row = open(tmp_path / "pd.csv").read().splitlines()[1].split(",")
    assert row[3:5] == ["5.0000", "6.0000"]
