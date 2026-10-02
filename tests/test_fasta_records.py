"""Each background-FASTA record is reported once, under a name of its own.

A FASTA can hold two records with one id (the version stripped), records whose header
gives no id, records with no sequence, and headers that are not GENCODE's or Ensembl's.
The background built from it must neither lose a competitor nor report one record twice:

- a record with a configured transcript's sequence and another with the same id went into
  two tables under one name, so the competitor was lost (when the id was the gene
  background's) or reported both as a copy and as a competitor;
- an empty record with a blank line after its header replaced a gene-background
  transcript with an empty sequence, which has no window: exit 2;
- three header-less copies were reported as one, and the competitor was named by its place
  among the records kept, not in the file;
- a header whose first ``|`` field is empty (``>|x``) raised IndexError, as in 2.4.1.
"""
import json
import random

import pytest

from isoform_dominance import cli
from isoform_dominance import identifiability as I


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


SHARED, ALT_A, ALT_B = _seq(900, 81), _seq(500, 82), _seq(600, 83)
T1, T2 = "ENST00000000001", "ENST00000000002"
SEQS = {T1: SHARED + ALT_A, T2: SHARED + ALT_B}
CFG = {"gene": "GENEX", "groups": {"A": [T1], "B": [T2]}, "primary_comparison": ["A", "B"]}
X = ALT_A[:250] + SHARED[300:700] + ALT_B[:200] + _seq(300, 84)   # competes with both
G = "ENST00000000007"                                              # a gene-background id
G_SEQ = _seq(800, 85)


def _run(tmp_path, capsys, fasta_text, background=None, *extra):
    cfg, sq, fa = tmp_path / "cfg.json", tmp_path / "seqs.json", tmp_path / "bg.fa"
    cfg.write_text(json.dumps(CFG))
    sq.write_text(json.dumps(SEQS))
    fa.write_text(fasta_text)
    argv = ["identifiability", "--config", str(cfg), "--sequences", str(sq),
            "--no-gene-background", "--background-fasta", str(fa), "--json", *extra]
    if background is not None:
        bs = tmp_path / "bg.json"
        bs.write_text(json.dumps(background))
        argv += ["--background-sequences", str(bs)]
    rc = cli.main(argv)
    out, err = capsys.readouterr()
    return rc, json.loads(out), err


def _answer(rc, res):
    return {"exit": rc, "verdict": res["verdict"], "rank": res["contrast"]["rank"],
            "unique": [res["groups"][g]["n_unique_kmers"] for g in ("A", "B")],
            "contrast": pytest.approx(res["contrast"]["min_resolvable_log2fc"], rel=1e-9),
            "A": pytest.approx(res["groups"]["A"]["min_resolvable_log2fc"], rel=1e-9)}


def _x_as_a_column(tmp_path, capsys):
    """The answer with X a background column and nothing else in the background."""
    rc, res, _ = _run(tmp_path, capsys, "", {"ENST00000000099": X})
    return _answer(rc, res)


@pytest.mark.parametrize("copy_first", [True, False])
def test_two_records_of_a_gene_id_one_a_configured_copy_lose_no_competitor(
        tmp_path, capsys, copy_first):
    recs = [">%s.1\n%s\n" % (G, SEQS[T1]), ">%s.2\n%s\n" % (G, X)]
    text = "".join(recs if copy_first else recs[::-1])
    rc, res, _ = _run(tmp_path, capsys, text, {G: G_SEQ})
    assert _answer(rc, res) == _x_as_a_column(tmp_path, capsys)
    bg = res["background"]
    names = set(bg["identical_to_configured"]) | set(bg["fasta_competitors"]) \
        | set(bg["gene_transcripts"])
    assert names == {G, G + "#2"}                    # two records, two names


@pytest.mark.parametrize("copy_first", [True, False])
def test_two_records_of_one_id_are_never_reported_both_ways(tmp_path, capsys, copy_first):
    y = "ENST00000000077"
    recs = [">%s\n%s\n" % (y, SEQS[T1]), ">%s\n%s\n" % (y, X)]
    rc, res, _ = _run(tmp_path, capsys, "".join(recs if copy_first else recs[::-1]))
    bg = res["background"]
    assert not set(bg["identical_to_configured"]) & set(bg["fasta_competitors"])
    assert len(bg["identical_to_configured"]) == len(bg["fasta_competitors"]) == 1
    assert _answer(rc, res) == _x_as_a_column(tmp_path, capsys)


@pytest.mark.parametrize("blank", ["\n", "\r\n", "  \n"])
def test_an_empty_record_with_a_blank_line_is_no_record(tmp_path, capsys, blank):
    rest = ">ENST00000000099\n%s\n" % X
    rc0, res0, _ = _run(tmp_path, capsys, ">%s\n%s" % (G, rest), {G: G_SEQ})
    rc, res, _ = _run(tmp_path, capsys, ">%s\n%s%s" % (G, blank, rest), {G: G_SEQ})
    assert res["background"]["sequence_from_fasta"] == []
    assert res["gene_total"]["estimable"] is True
    assert _answer(rc, res) == _answer(rc0, res0)


def test_header_less_records_are_each_named_by_their_place_in_the_file(tmp_path, capsys):
    text = (">\n%s\n>\n%s\n>\n%s\n>ENST00000000098\n%s\n>\n%s\n"
            % (SEQS[T1], SEQS[T2], SEQS[T1], _seq(400, 86), X))
    rc, res, _ = _run(tmp_path, capsys, text)
    bg = res["background"]
    assert bg["identical_to_configured"] == {"record1": T1, "record2": T2, "record3": T1}
    assert list(bg["fasta_competitors"]) == ["record5"]


@pytest.mark.parametrize("head", ["|x", " |x", "|", ".1", ">"])
def test_a_header_with_an_empty_first_field_is_read(tmp_path, head):
    path = tmp_path / "bg.fa"
    path.write_text(">%s\n%s\n" % (head, X))
    query = I.kmers(SEQS[T1], 31)
    assert I.scan_background_fasta(path, query, 31) == I.kmers(X, 31) & query
    got = I.scan_fasta_competitors(path, query, 31)
    assert len(got) == 1
