"""A decoy-aware index's gentrome as --background-fasta.

The FASTA a decoy-aware Salmon index is built from holds the genome after the transcripts.
Each exon is in the genome, so every window inside one has a copy there and, judged against
it, loses its uniqueness, and since 2.5.0 the genome record would be a column of the system
too.  Salmon sets aside only the reads that map better to a decoy than to any transcript;
an exon's reads match the transcript as well as the genome, so a decoy competes for none of
them.  ``--decoys decoys.txt`` skips those records; without it a record longer than 1 Mb is
named in a NOTE.
"""
import json
import random

import pytest

from isoform_dominance import cli
from isoform_dominance import identifiability as I


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


E1, E2, E3, E4 = _seq(600, 61), _seq(300, 62), _seq(350, 63), _seq(900, 64)
T1, T2 = "ENST00000000001", "ENST00000000002"
SEQS = {T1: E1 + E2 + E4, T2: E1 + E3 + E4}
OTHER = _seq(1200, 65)
# the gene's exons in the genome, with introns and intergenic sequence around them
CHROM = (_seq(9000, 66) + E1 + _seq(2000, 67) + E2 + _seq(1500, 68) + E3 + _seq(2500, 69)
         + E4 + _seq(12000, 70))
CFG = {"gene": "GENEX", "groups": {"A": [T1], "B": [T2]}, "primary_comparison": ["A", "B"]}


def _gencode(tid, seq):
    return ">%s.1|ENSG00000000010.3|-|-|GENEX-201|GENEX|%d|protein_coding|\n%s\n" % (
        tid, len(seq), seq)


TRANSCRIPTS = (_gencode(T1, SEQS[T1]) + _gencode(T2, SEQS[T2])
               + _gencode("ENST00000000009", OTHER))
GENTROME = TRANSCRIPTS + ">chrF\n%s\n>chrG extra words\n%s\n" % (CHROM, _seq(500, 71))


@pytest.fixture
def small_genome(monkeypatch):
    """The chromosome above is 30 kb; a record longer than 20 kb counts as long here."""
    monkeypatch.setattr(I, "LONG_RECORD", 20_000)


def _run(tmp_path, capsys, fasta, *extra):
    cfg, sq, fa = tmp_path / "cfg.json", tmp_path / "seqs.json", tmp_path / "bg.fa"
    cfg.write_text(json.dumps(CFG))
    sq.write_text(json.dumps(SEQS))
    fa.write_text(fasta)
    rc = cli.main(["identifiability", "--config", str(cfg), "--sequences", str(sq),
                   "--no-gene-background", "--background-fasta", str(fa), "--min-log2fc",
                   "0.5", "--json", *extra])
    out, err = capsys.readouterr()
    return rc, json.loads(out), err


def _decoys(tmp_path, *names):
    p = tmp_path / "decoys.txt"
    p.write_text("".join(n + "\n" for n in names))
    return str(p)


def _answer(rc, res):
    return {"exit": rc, "verdict": res["verdict"], "rank": res["contrast"]["rank"],
            "contrast": res["contrast"]["min_resolvable_log2fc"],
            "unique": {g: res["groups"][g]["n_unique_kmers"] for g in ("A", "B")}}


def test_the_threshold_is_1_mb():
    assert I.LONG_RECORD == 1_000_000


def test_with_decoys_the_gentrome_gives_the_transcripts_answer(tmp_path, capsys, small_genome):
    alone = _run(tmp_path, capsys, TRANSCRIPTS)
    rc, res, err = _run(tmp_path, capsys, GENTROME, "--decoys",
                        _decoys(tmp_path, "chrF", "chrG"))
    assert _answer(rc, res) == _answer(*alone[:2])
    bg = res["background"]
    assert (bg["decoys_listed"], bg["decoys_skipped"], bg["decoys_absent"]) == (2, 2, [])
    assert bg["fasta_long_records"] == {} and bg["fasta_competitors"] == {}
    assert "genome" not in err


def test_without_decoys_the_genome_is_named_and_competes(tmp_path, capsys, small_genome):
    alone = _run(tmp_path, capsys, TRANSCRIPTS)
    rc, res, err = _run(tmp_path, capsys, GENTROME)
    bg = res["background"]
    assert bg["fasta_long_records"] == {"chrF": len(CHROM)}
    assert "chrF" in bg["fasta_competitors"]
    # every window inside an exon now has a copy: only the junctions stay unique
    assert 0 < res["groups"]["A"]["n_unique_kmers"] < alone[1]["groups"]["A"]["n_unique_kmers"]
    assert "looks like genome sequence" in err and "--decoys decoys.txt" in err
    assert "without the genome decoys" in err


def test_a_long_record_decoys_does_not_name_is_said_too(tmp_path, capsys, small_genome):
    rc, res, err = _run(tmp_path, capsys, GENTROME, "--decoys", _decoys(tmp_path, "chrG"))
    assert res["background"]["fasta_long_records"] == {"chrF": len(CHROM)}
    assert "--decoys %s does not name (chrF)" % (tmp_path / "decoys.txt") in err


def test_decoy_names_the_fasta_does_not_hold_are_named(tmp_path, capsys, small_genome):
    rc, res, err = _run(tmp_path, capsys, GENTROME, "--decoys",
                        _decoys(tmp_path, "chrF", "chrG", "chr21"))
    assert res["background"]["decoys_absent"] == ["chr21"]
    assert "does not hold (chr21)" in err


def test_decoys_without_a_background_fasta_is_refused(tmp_path, capsys):
    cfg, sq = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps(CFG))
    sq.write_text(json.dumps(SEQS))
    rc = cli.main(["identifiability", "--config", str(cfg), "--sequences", str(sq),
                   "--decoys", _decoys(tmp_path, "chrF")])
    assert rc == 1 and "--decoys" in capsys.readouterr().err


def test_the_harness_call_takes_decoys(tmp_path, small_genome):
    fa = tmp_path / "bg.fa"
    fa.write_text(GENTROME)
    kw = {"k": 31, "window": 31, "sequences": SEQS, "background_gene_transcripts": False,
          "min_log2fc": 0.5}
    fa2 = tmp_path / "tx.fa"
    fa2.write_text(TRANSCRIPTS)
    with_decoys = I.analyze(CFG, background_fasta=str(fa),
                            decoys=_decoys(tmp_path, "chrF", "chrG"), **kw)
    alone = I.analyze(CFG, background_fasta=str(fa2), **kw)
    assert with_decoys["contrast"]["min_resolvable_log2fc"] \
        == alone["contrast"]["min_resolvable_log2fc"]


def test_a_decoy_is_matched_by_its_first_word_or_that_words_first_field(tmp_path):
    fa = tmp_path / "bg.fa"
    fa.write_text(">chrF|whole chromosome\n%s\n>chrG\n%s\n" % (CHROM, E2))
    stats = {}
    got = I.scan_fasta_competitors(str(fa), I.kmers(SEQS[T1], 31), 31,
                                   decoys=["chrF"], stats=stats)
    assert list(got) == ["chrG"] and stats["decoys_skipped"] == 1
    assert stats["decoys_found"] == {"chrF"}
