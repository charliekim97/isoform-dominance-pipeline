"""Identical sequences in the background are one sequence, wherever they come from.

Salmon's default index keeps one of identical sequences.  2.4.1 applied that to a
``--background-fasta`` record with a configured transcript's sequence, and to nothing else:
a transcript of the gene background with a configured transcript's sequence still took
every window that transcript had (NR1H3 at release 110 is ``not_identifiable`` for that
reason alone), and a FASTA record with a gene-background transcript's sequence -- the chrY
copy of a pseudoautosomal transcript -- was warned about as a copy that splits reads.  The
rule now applies to the whole background: the gene's, ``background_sequences``,
``--inputs`` and the FASTA's records.  ``--keep-duplicates`` counts every copy, as before.
"""
import json
import random
import urllib.request

import pytest

from isoform_dominance import cli
from isoform_dominance import identifiability as I

SERVER = "https://rest.ensembl.org"


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


SHARED, ALT_A, ALT_B, OWN = _seq(1500, 51), _seq(800, 52), _seq(1000, 53), _seq(700, 54)
T1, T2 = "ENST00000000001", "ENST00000000002"
SEQS = {T1: SHARED + ALT_A, T2: SHARED + ALT_B}
TWIN = "ENST00000000003"                  # a gene-background transcript with T1's sequence
OTHER = "ENST00000000004"                 # a gene-background transcript of its own
BG = {TWIN: SEQS[T1], OTHER: SHARED[:900] + OWN}
GENE_ID = "ENSG00000000010"
CFG = {"gene": "GENEX", "gene_id": GENE_ID, "groups": {"A": [T1], "B": [T2]},
       "primary_comparison": ["A", "B"]}


def _answer(res, rc=None):
    out = {"rank": res["contrast"]["rank"], "classes": res["n_compatibility_classes"],
           "contrast": res["contrast"]["min_resolvable_log2fc"],
           "verdict": res["verdict"], "effect": res["effect_resolvable"]}
    for g in ("A", "B"):
        out[g] = (res["groups"][g]["n_unique_kmers"], res["groups"][g]["min_resolvable_log2fc"])
    if rc is not None:
        out["exit"] = rc
    return out


# --------------------------------------------------------------------------- #
# a gene-background twin of a configured transcript, from each source
# --------------------------------------------------------------------------- #
def test_background_sequences_twin_is_not_a_competitor():
    kw = {"sequences": SEQS, "background_gene_transcripts": False, "min_log2fc": 0.5}
    with_twin = I.analyze(CFG, background_sequences=BG, **kw)
    without = I.analyze(CFG, background_sequences={OTHER: BG[OTHER]}, **kw)
    assert _answer(with_twin) == _answer(without)
    assert with_twin["groups"]["A"]["n_unique_kmers"] > 0
    bg = with_twin["background"]
    assert bg["identical_to_configured"] == {TWIN: T1}
    assert bg["identical_source"] == {TWIN: "sequences"}
    assert bg["gene_transcripts"] == [OTHER] and bg["n_background_transcripts"] == 1


def test_keep_duplicates_counts_the_twin_as_2_4_1_did():
    res = I.analyze(CFG, sequences=SEQS, background_sequences=BG,
                    background_gene_transcripts=False, keep_duplicates=True)
    assert res["groups"]["A"]["n_unique_kmers"] == 0
    assert res["background"]["identical_to_configured"] is None
    assert res["background"]["gene_transcripts"] == sorted([OTHER, TWIN])


def _files(tmp_path, bg=None):
    cfg, sq = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps(CFG))
    sq.write_text(json.dumps(SEQS))
    argv = ["identifiability", "--config", str(cfg), "--sequences", str(sq),
            "--min-log2fc", "0.5", "--json"]
    if bg is not None:
        bs = tmp_path / "bg.json"
        bs.write_text(json.dumps(bg))
        argv += ["--background-sequences", str(bs)]
    return argv


def _run(argv, capsys):
    rc = cli.main(argv)
    out, err = capsys.readouterr()
    return rc, json.loads(out), err


def test_cli_background_sequences_twin_is_a_note_not_a_competitor(tmp_path, capsys):
    rc, res, err = _run(_files(tmp_path, BG), capsys)
    rc0, res0, _ = _run(_files(tmp_path, {OTHER: BG[OTHER]}), capsys)
    assert _answer(res, rc) == _answer(res0, rc0)
    assert "NOTE: the gene background: 1 transcript(s) with the sequence of a configured " \
           "transcript were not counted" in err
    assert "%s (= %s)" % (TWIN, T1) in err
    # and what that means for extract, with the two ways out
    assert "comes first" in err and "`extract`" in err
    assert "same group" in err and "--keepDuplicates" in err


class _Rest:
    """The gene as Ensembl gives it: the configured transcripts and ``background``."""

    def __init__(self, background):
        self.seqs = dict(SEQS, **background)

    def __call__(self, req, timeout=None, **_):
        from test_ensembl_http import _Resp
        path = req.full_url[len(SERVER):]
        if path.startswith("/info/data"):
            return _Resp(json.dumps({"releases": [116]}))
        if path.startswith("/lookup/id/%s" % GENE_ID):
            return _Resp(json.dumps({"id": GENE_ID, "Transcript": [
                {"id": t + ".1"} for t in sorted(self.seqs)]}))
        if path.startswith("/sequence/id") and req.get_method() == "POST":
            return _Resp(json.dumps([{"query": i, "id": i, "seq": self.seqs[i.split(".")[0]]}
                                     for i in json.loads(req.data)["ids"]]))
        raise AssertionError("unexpected request %s" % path)


def _fetched(monkeypatch, background, **kw):
    monkeypatch.setattr(urllib.request, "urlopen", _Rest(background))
    return I.analyze(CFG, min_log2fc=0.5, **kw)


def test_a_fetched_twin_is_not_a_competitor(monkeypatch):
    with_twin = _fetched(monkeypatch, BG)
    without = _fetched(monkeypatch, {OTHER: BG[OTHER]})
    assert _answer(with_twin) == _answer(without)
    assert with_twin["background"]["identical_to_configured"] == {TWIN: T1}
    assert with_twin["background"]["identical_source"] == {TWIN: "gene"}
    counted = _fetched(monkeypatch, BG, keep_duplicates=True)
    assert counted["groups"]["A"]["n_unique_kmers"] == 0


def test_a_rerun_from_inputs_leaves_the_twin_out_and_says_where_it_came_from(
        monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(urllib.request, "urlopen", _Rest(BG))
    cfg, saved = tmp_path / "cfg.json", tmp_path / "in.json"
    cfg.write_text(json.dumps(CFG))
    rc, live, _ = _run(["identifiability", "--config", str(cfg), "--min-log2fc", "0.5",
                        "--save-inputs", str(saved), "--json"], capsys)
    assert json.loads(saved.read_text())["background_sequences"][TWIN] == SEQS[T1]

    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    rc2, rerun, _ = _run(["identifiability", "--config", str(cfg), "--inputs", str(saved),
                          "--min-log2fc", "0.5", "--json"], capsys)
    assert (rc2, _answer(rerun)) == (rc, _answer(live))
    for key in ("identical_to_configured", "identical_source", "gene_transcripts"):
        assert rerun["background"][key] == live["background"][key], key
    assert live["background"]["identical_source"] == {TWIN: "gene"}


def test_saved_inputs_record_keep_duplicates_for_a_gene_background(
        monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(urllib.request, "urlopen", _Rest(BG))
    cfg, saved = tmp_path / "cfg.json", tmp_path / "in.json"
    cfg.write_text(json.dumps(CFG))
    cli.main(["identifiability", "--config", str(cfg), "--keep-duplicates",
              "--save-inputs", str(saved)])
    assert json.loads(saved.read_text())["analysis"]["keep_duplicates"] is True
    capsys.readouterr()
    cli.main(["identifiability", "--config", str(cfg), "--inputs", str(saved)])
    assert ("the inputs were saved with --keep-duplicates and this run is without it"
            in capsys.readouterr().err)


# --------------------------------------------------------------------------- #
# a FASTA record with a gene-background transcript's sequence (B5b)
# --------------------------------------------------------------------------- #
def _gencode(tid, gid, seq):
    return ">%s.1|%s.3|-|-|GENEX-201|GENEX|%d|protein_coding|\n%s\n" % (tid, gid, len(seq), seq)


CHRY = "ENST00000000066"                  # the chrY copy of OTHER: its own gene id
CHRY_GENE = "ENSG00000000060"
REF = "".join(_gencode(t, GENE_ID, s) for t, s in sorted(dict(SEQS, **{OTHER: BG[OTHER]})
                                                         .items()))


def _with_fasta(tmp_path, fasta_text, *extra):
    fa = tmp_path / "bg.fa"
    fa.write_text(fasta_text)
    return _files(tmp_path, {OTHER: BG[OTHER]}) + ["--background-fasta", str(fa), *extra]


def test_a_fasta_copy_of_a_background_transcript_is_that_transcript(tmp_path, capsys):
    rc, res, err = _run(_with_fasta(tmp_path, REF + _gencode(CHRY, CHRY_GENE, BG[OTHER])),
                        capsys)
    rc0, res0, _ = _run(_with_fasta(tmp_path, REF), capsys)
    assert _answer(res, rc) == _answer(res0, rc0)
    bg = res["background"]
    assert bg["identical_to_background"] == {CHRY: OTHER}
    assert bg["identical_source"] == {CHRY: "fasta"}
    assert bg["fasta_competitors"] == {}
    assert [c["gene_id"] for c in bg["same_name_copies"]] == [CHRY_GENE]   # still listed
    assert "WARNING" not in err and "reads are split" not in err
    assert "%s (= %s)" % (CHRY, OTHER) in err and "counted once" in err


def test_a_copy_one_base_different_is_warned_about_and_competes(tmp_path, capsys):
    s = BG[OTHER]
    off = s[:300] + ("A" if s[300] != "A" else "C") + s[301:]
    rc, res, err = _run(_with_fasta(tmp_path, REF + _gencode(CHRY, CHRY_GENE, off)), capsys)
    assert res["background"]["identical_to_background"] == {}
    assert list(res["background"]["fasta_competitors"]) == [CHRY]
    assert "WARNING" in err and "reads are split" in err


def test_keep_duplicates_warns_and_counts_the_copy(tmp_path, capsys):
    rc, res, err = _run(_with_fasta(tmp_path, REF + _gencode(CHRY, CHRY_GENE, BG[OTHER]),
                                    "--keep-duplicates"), capsys)
    assert res["background"]["identical_to_background"] is None
    assert list(res["background"]["fasta_competitors"]) == [CHRY]
    assert "WARNING" in err and "reads are split" in err


@pytest.mark.parametrize("keep", [False, True])
def test_identical_fasta_records_are_one_column(tmp_path, capsys, keep):
    r1, r2 = "ENST00000000071", "ENST00000000072"
    copy = ALT_B[:500] + OWN
    text = REF + _gencode(r1, "ENSG00000000070", copy) + _gencode(r2, "ENSG00000000070", copy)
    rc, res, err = _run(_with_fasta(tmp_path, text, *(["--keep-duplicates"] if keep else [])),
                        capsys)
    bg = res["background"]
    if keep:
        assert sorted(bg["fasta_competitors"]) == [r1, r2]
    else:
        assert list(bg["fasta_competitors"]) == [r1]
        assert bg["identical_to_background"] == {r2: r1}


def test_a_copy_of_a_background_transcript_that_shares_nothing_is_not_warned_about(
        tmp_path, capsys):
    """The scan keeps a record that shares no window only to tell that it is a copy: it is
    no column, and with a gene-background transcript's sequence no reason to warn either."""
    apart = "ENST00000000005"
    alone = _seq(600, 55)                      # shares no window with either class
    bg = {OTHER: BG[OTHER], apart: alone}
    fa = tmp_path / "bg.fa"
    fa.write_text(REF + _gencode(apart, GENE_ID, alone) + _gencode(CHRY, CHRY_GENE, alone))
    rc, res, err = _run(_files(tmp_path, bg) + ["--background-fasta", str(fa)], capsys)
    assert res["background"]["identical_to_background"] == {CHRY: apart}
    assert res["background"]["fasta_competitors"] == {}
    assert "WARNING" not in err


def test_a_gene_transcript_the_fasta_holds_as_a_configured_ones_sequence(tmp_path, capsys):
    """One id in both backgrounds, and the FASTA's sequence for it is a configured
    transcript's: the FASTA's is the one the index has, so it is a copy of that transcript,
    not the gene background's competitor."""
    fa = tmp_path / "bg.fa"
    fa.write_text(_gencode(OTHER, GENE_ID, SEQS[T1]))
    rc, res, err = _run(_files(tmp_path, {OTHER: BG[OTHER]}) + ["--background-fasta", str(fa)],
                        capsys)
    bg = res["background"]
    assert bg["sequence_from_fasta"] == [OTHER]
    assert bg["identical_to_configured"] == {OTHER: T1}
    assert bg["identical_source"] == {OTHER: "fasta"}
    assert bg["gene_transcripts"] == []
    rc0, res0, _ = _run(_files(tmp_path), capsys)          # no background at all
    assert _answer(res, rc) == _answer(res0, rc0)


def test_a_copy_whose_header_gives_no_id_is_still_reported(tmp_path, capsys):
    fa = tmp_path / "bg.fa"
    fa.write_text(">\n%s\n>%s\n%s\n" % (SEQS[T1], CHRY, SEQS[T1]))
    rc, res, err = _run(_files(tmp_path, {OTHER: BG[OTHER]}) + ["--background-fasta", str(fa)],
                        capsys)
    assert res["background"]["identical_to_configured"] == {"(no id 1)": T1, CHRY: T1}
