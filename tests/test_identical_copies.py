"""--background-fasta and sequences identical to a configured transcript.

Salmon's index keeps only the first of identical sequences unless it was built with
``--keepDuplicates``: in GENCODE v50's reference-chromosome index it removed 1,600, among
them all 382 chrY transcripts of the 18 protein-coding genes annotated on both chrX and
chrY.  The FASTA an index is built from still holds them, and a background FASTA record
with a configured transcript's sequence took every unique k-mer that transcript had: a
CD99 config of the chrX gene had none in either class against a FASTA of CD99's chrX and
chrY records, and 1,273 and 3,051 without it.
"""
import json
import random

from isoform_dominance import cli, io


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


_SHARED = _seq(1500, 41)
_ALT_A = _seq(800, 42)
_ALT_B = _seq(1000, 43)
T1, T2 = "ENST00000000001", "ENST00000000002"
SEQS = {T1: _SHARED + _ALT_A, T2: _SHARED + _ALT_B}
GROUPS = {"A": [T1], "B": [T2]}

TWIN = "ENST00000000009"                         # T1's sequence under another id


def _one_off(seq, at):
    return seq[:at] + {"A": "C", "C": "G", "G": "T", "T": "A"}[seq[at]] + seq[at + 1:]


def _case(tmp_path, fasta_records=None, groups=GROUPS, seqs=SEQS, name="bg.fa"):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"gene": "GENEX", "groups": groups,
                               "primary_comparison": ["A", "B"]}))
    sq = tmp_path / "seqs.json"
    sq.write_text(json.dumps(seqs))
    argv = ["identifiability", "--config", str(cfg), "--sequences", str(sq),
            "--no-gene-background", "--json"]
    if fasta_records is not None:
        fa = tmp_path / name
        fa.write_text("".join(">%s|ENSG00000000077.1|-|-|X-201|X|%d|protein_coding|\n%s\n"
                              % (t, len(s), s) for t, s in fasta_records.items()))
        argv += ["--background-fasta", str(fa)]
    return argv


def _run(argv, capsys):
    rc = cli.main(argv)
    out, err = capsys.readouterr()
    return rc, json.loads(out), err


def _unique(res):
    return {g: e["n_unique_kmers"] for g, e in res["groups"].items()}


def test_an_identical_record_is_not_counted_by_default(tmp_path, capsys):
    _, alone, _ = _run(_case(tmp_path), capsys)
    rc, res, err = _run(_case(tmp_path, {TWIN: SEQS[T1]}), capsys)
    assert _unique(alone)["A"] > 0
    assert _unique(res) == _unique(alone)
    assert res["verdict"] == alone["verdict"] and rc == cli.EXIT_OK
    assert res["background"]["identical_to_configured"] == {TWIN: T1}
    assert res["background"]["keep_duplicates"] is False
    assert "Salmon's default index keeps one of identical sequences" in err
    assert TWIN in err and "--keep-duplicates" in err


def test_keep_duplicates_counts_it_as_before(tmp_path, capsys):
    rc, res, err = _run(_case(tmp_path, {TWIN: SEQS[T1]}) + ["--keep-duplicates"], capsys)
    assert _unique(res)["A"] == 0
    assert _unique(res)["B"] > 0
    assert res["background"]["identical_to_configured"] is None      # not looked for
    assert res["background"]["keep_duplicates"] is True
    assert "identical sequences" not in err


def test_one_base_different_is_not_left_out(tmp_path, capsys):
    _, alone, _ = _run(_case(tmp_path), capsys)
    near = _one_off(SEQS[T1], len(_SHARED) + 400)
    _, res, err = _run(_case(tmp_path, {TWIN: near}), capsys)
    assert res["background"]["identical_to_configured"] == {}
    assert 0 < _unique(res)["A"] < _unique(alone)["A"]     # only windows over that base
    assert "identical sequences" not in err


def test_a_difference_of_case_or_trailing_space_is_no_difference(tmp_path, capsys):
    _, alone, _ = _run(_case(tmp_path), capsys)
    fa = tmp_path / "bg.fa"
    fa.write_text(">%s|ENSG00000000077.1|-|-|X-201|X|%d|protein_coding|\n%s   \n"
                  % (TWIN, len(SEQS[T1]), SEQS[T1].lower()))
    argv = _case(tmp_path) + ["--background-fasta", str(fa)]
    _, res, _ = _run(argv, capsys)
    assert res["background"]["identical_to_configured"] == {TWIN: T1}
    assert _unique(res) == _unique(alone)


def test_two_configured_transcripts_with_one_sequence_stay_indistinguishable(
        tmp_path, capsys):
    """Only background records are left out.  Two configured transcripts that are one
    sequence, in different classes, cannot be told apart -- with or without a FASTA."""
    seqs = {T1: SEQS[T1], T2: SEQS[T1]}
    for fasta in (None, {TWIN: SEQS[T1]}):
        _, res, _ = _run(_case(tmp_path, fasta, seqs=seqs), capsys)
        assert _unique(res) == {"A": 0, "B": 0}


def test_only_a_record_with_a_configured_sequence_is_left_out(tmp_path, capsys):
    """A background record identical to another background record, or to nothing, stays."""
    other = _seq(900, 44)
    _, res, _ = _run(_case(tmp_path, {TWIN: SEQS[T1], "ENST00000000010": other,
                                      "ENST00000000011": other}), capsys)
    assert res["background"]["identical_to_configured"] == {TWIN: T1}


def test_saved_inputs_record_the_setting_and_a_rerun_with_another_says_so(
        tmp_path, capsys):
    base = _case(tmp_path, {TWIN: SEQS[T1]})
    saved = tmp_path / "in.json"
    for flag, expect in (([], False), (["--keep-duplicates"], True)):
        cli.main(base + flag + ["--save-inputs", str(saved)])
        capsys.readouterr()
        assert io.load_inputs(str(saved))["analysis"]["keep_duplicates"] is expect
    # saved with --keep-duplicates; rerun without it
    cfg = base[base.index("--config") + 1]
    fa = base[base.index("--background-fasta") + 1]
    cli.main(["identifiability", "--config", cfg, "--inputs", str(saved),
              "--background-fasta", fa, "--json"])
    err = capsys.readouterr().err
    # the note on the rerun itself, not the one about left-out records, which names the
    # flag too
    assert "the inputs were saved with --keep-duplicates and this run is without it" in err
