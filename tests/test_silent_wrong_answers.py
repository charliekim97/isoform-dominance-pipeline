"""Two configurations that gave an answer, exit 0, and the answer was wrong.

A window longer than the read: no read can hold a whole window, so no fragment is ever
informative.  The informative fraction fell from 0.119 at window 100 to 0.0 at 101 with
100-nt reads, and the run reported that as a property of the gene.

A transcript in two classes: its +1 and -1 in the contrast cancel when the two are the
compared pair, and in any two its column entered the system twice with the first copy all
zero.  The CLI then reported a precondition failure -- "A transcript is shorter than
window=31", naming none -- and exit 2, for transcripts thousands of bases long.
"""
import json
import random

import pytest

from isoform_dominance import cli, identifiability


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


_SHARED = _seq(1500, 41)
SEQS = {"A1": _SHARED + _seq(600, 42), "B1": _SHARED + _seq(700, 43), "C1": _seq(900, 44)}


def _case(tmp_path, groups, primary=("A", "B")):
    cfg, sq = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps({"groups": groups, "primary_comparison": list(primary)}))
    sq.write_text(json.dumps(SEQS))
    return ["identifiability", "--config", str(cfg), "--sequences", str(sq),
            "--no-gene-background"]


def _err_line(capsys):
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1, err
    return err[0]


@pytest.mark.parametrize("window,reads,ok", [(100, 100, True), (101, 100, False),
                                             (76, 75, False), (75, 75, True)])
def test_a_window_longer_than_the_read_is_refused(tmp_path, capsys, window, reads, ok):
    argv = _case(tmp_path, {"A": ["A1"], "B": ["B1"]}) + [
        "--window", str(window), "--read-length", str(reads)]
    code = cli.main(argv)
    if ok:
        assert code == cli.EXIT_OK
        capsys.readouterr()
    else:
        assert code == 1
        assert "window %d exceeds the read length %d" % (window, reads) in _err_line(capsys)


def test_the_library_refuses_it_too():
    with pytest.raises(ValueError, match="exceeds the read length"):
        identifiability.analyze({"groups": {"A": ["A1"], "B": ["B1"]}},
                                sequences=SEQS, background_gene_transcripts=False,
                                window=151, read_length=150)


def test_a_transcript_in_both_compared_classes_is_refused(tmp_path, capsys):
    argv = _case(tmp_path, {"A": ["A1", "C1"], "B": ["B1", "C1"]})
    assert cli.main(argv) == 1
    line = _err_line(capsys)
    assert "C1" in line and '"A"' in line and '"B"' in line


def test_a_transcript_shared_with_a_third_class_is_refused_too(tmp_path, capsys):
    # not a cancelling contrast, but the same duplicated column
    argv = _case(tmp_path, {"A": ["A1"], "B": ["B1", "C1"], "C": ["C1"]})
    assert cli.main(argv) == 1
    line = _err_line(capsys)
    assert "C1" in line and '"B"' in line and '"C"' in line


def test_a_versioned_and_an_unversioned_id_are_one_transcript(tmp_path, capsys):
    argv = _case(tmp_path, {"A": ["A1", "C1.2"], "B": ["B1", "C1"]})
    assert cli.main(argv) == 1
    assert "C1" in _err_line(capsys)
