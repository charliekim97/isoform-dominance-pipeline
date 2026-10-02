"""A competing sequence gives one answer whichever way it is passed.

Through 2.4.1 a ``--background-fasta`` record took windows from the uniqueness layer only:
the compatibility system was built from the gene background alone.  The same competitor
passed with ``--background-sequences`` was a column of that system.  So the recommended
route gave the optimistic answer: one fewer rank, a smaller min |log2FC|, and exit 0 where
the other route exits 3.  A record that shares a window with a configured transcript is
now a column too.
"""
import json
import random

import pytest

from isoform_dominance import cli
from isoform_dominance import identifiability as I


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


SHARED, ALT_A, ALT_B, OWN = _seq(1500, 1), _seq(400, 2), _seq(400, 3), _seq(300, 4)
T1, T2, COMP = "ENST00000000001", "ENST00000000002", "ENST00000000099"
SEQS = {T1: SHARED + ALT_A, T2: SHARED + ALT_B}
# most of A's own exon, a stretch of the shared one, and sequence of its own
COMPETITOR = ALT_A[:330] + SHARED[:600] + OWN
CFG = {"gene": "GENEX", "groups": {"A": [T1], "B": [T2]}, "primary_comparison": ["A", "B"]}
STRANGER = "ENST00000000098"                  # shares no window with either class


def _fasta(path, records):
    path.write_text("".join(">%s.1|ENSG00000000077.1|-|-|Y-201|Y|%d|protein_coding|\n%s\n"
                            % (t, len(s), s) for t, s in records.items()))
    return str(path)


def _cli(tmp_path, capsys, *extra):
    cfg, sq = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps(CFG))
    sq.write_text(json.dumps(SEQS))
    rc = cli.main(["identifiability", "--config", str(cfg), "--sequences", str(sq),
                   "--no-gene-background", "--min-log2fc", "0.5", "--json", *extra])
    out, err = capsys.readouterr()
    return rc, json.loads(out), err


def _answer(rc, res):
    return {"rank": res["contrast"]["rank"],
            "contrast": res["contrast"]["min_resolvable_log2fc"],
            "A": res["groups"]["A"]["min_resolvable_log2fc"],
            "B": res["groups"]["B"]["min_resolvable_log2fc"],
            "classes": res["n_compatibility_classes"], "exit": rc}


def test_fasta_and_background_sequences_give_one_answer(tmp_path, capsys):
    bs = tmp_path / "bg.json"
    bs.write_text(json.dumps({COMP: COMPETITOR}))
    fa = _fasta(tmp_path / "bg.fa", {COMP: COMPETITOR})
    by_seqs = _answer(*_cli(tmp_path, capsys, "--background-sequences", str(bs))[:2])
    by_fasta = _answer(*_cli(tmp_path, capsys, "--background-fasta", fa)[:2])
    alone = _answer(*_cli(tmp_path, capsys)[:2])
    assert by_fasta == by_seqs
    # and the competitor matters, so the equality is not two runs that both ignore it
    assert by_seqs["rank"] == 3 and alone["rank"] == 2
    assert by_seqs["exit"] == cli.EXIT_EFFECT_NOT_RESOLVED and alone["exit"] == cli.EXIT_OK
    assert by_seqs["contrast"] > 0.5 > alone["contrast"]


def test_the_harness_call_gives_one_answer_too(tmp_path):
    fa = _fasta(tmp_path / "bg.fa", {COMP: COMPETITOR})
    kw = {"k": 31, "window": 31, "sequences": SEQS, "background_gene_transcripts": False,
          "min_log2fc": 0.5}
    by_seqs = I.analyze(CFG, background_sequences={COMP: COMPETITOR}, **kw)
    by_fasta = I.analyze(CFG, background_fasta=fa, **kw)
    assert _answer(cli._identifiability_exit(by_fasta), by_fasta) \
        == _answer(cli._identifiability_exit(by_seqs), by_seqs)


def test_the_report_names_the_competitors_and_counts_them_apart(tmp_path, capsys):
    fa = _fasta(tmp_path / "bg.fa", {COMP: COMPETITOR, STRANGER: _seq(900, 5)})
    rc, res, err = _cli(tmp_path, capsys, "--background-fasta", fa)
    bg = res["background"]
    shared = I.kmers(COMPETITOR, 31) & (I.kmers(SEQS[T1], 31) | I.kmers(SEQS[T2], 31))
    assert bg["fasta_competitors"] == {COMP: len(shared)}       # distinct windows
    assert bg["n_fasta_competitors"] == 1
    assert bg["n_background_transcripts"] == 0                  # the gene's: none here
    assert bg["gene_transcripts"] == []


def test_a_record_that_shares_no_window_changes_nothing(tmp_path, capsys):
    fa1 = _fasta(tmp_path / "one.fa", {COMP: COMPETITOR})
    fa2 = _fasta(tmp_path / "two.fa", {COMP: COMPETITOR, STRANGER: _seq(900, 5)})
    one = _cli(tmp_path, capsys, "--background-fasta", fa1)
    two = _cli(tmp_path, capsys, "--background-fasta", fa2)
    assert _answer(*one[:2]) == _answer(*two[:2])
    assert STRANGER not in two[1]["background"]["fasta_competitors"]


@pytest.mark.parametrize("differ", [False, True])
def test_an_id_in_both_backgrounds_is_one_column(tmp_path, capsys, differ):
    """One transcript in the gene background and in the FASTA is counted once; when the two
    carry different sequence (another release), the FASTA's is used, with a NOTE."""
    other = ALT_A[:200] + SHARED[:500] + OWN          # what another release holds
    bs = tmp_path / "bg.json"
    bs.write_text(json.dumps({COMP + ".3": COMPETITOR}))
    fa = _fasta(tmp_path / "bg.fa", {COMP: other if differ else COMPETITOR})
    rc, res, err = _cli(tmp_path, capsys, "--background-sequences", str(bs),
                        "--background-fasta", fa)
    bg = res["background"]
    assert bg["gene_transcripts"] == [COMP] and bg["fasta_competitors"] == {}
    assert bg["sequence_from_fasta"] == ([COMP] if differ else [])
    # the answer is the one for the sequence that was used
    used = tmp_path / "used.json"
    used.write_text(json.dumps({COMP: other if differ else COMPETITOR}))
    alone = _cli(tmp_path, capsys, "--background-sequences", str(used))
    assert _answer(rc, res) == _answer(*alone[:2])
    assert ("taken from --background-fasta" in err) is differ


def test_a_system_too_large_for_memory_is_one_line_and_exit_1(tmp_path, capsys, monkeypatch):
    """Every record that shares a window is a column, so a genome record or a repeat many
    records share can make a system numpy cannot allocate: say so, not a traceback."""
    def too_large(*a, **k):
        raise MemoryError("Unable to allocate 3.73 GiB for an array with shape (8371, 59840)")
    monkeypatch.setattr(I, "compatibility_matrix", too_large)
    fa = _fasta(tmp_path / "bg.fa", {COMP: COMPETITOR})
    cfg, sq = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps(CFG))
    sq.write_text(json.dumps(SEQS))
    rc = cli.main(["identifiability", "--config", str(cfg), "--sequences", str(sq),
                   "--background-fasta", fa])
    err = capsys.readouterr().err
    assert rc == 1 and len(err.strip().splitlines()) == 1
    assert err.startswith("error: out of memory (Unable to allocate 3.73 GiB")
    assert "--max-window-records" in err and "--decoys" in err
