"""Without --background-fasta the NOTE gives the figure the 2.5.0 decision rule measured (issue #3).

The rule posted on issue #3 before the code was written: if the whole-index background changed
the exit status at --min-log2fc 0.5 for 10 or more of the 100 two-class configs annotate proposes
at Ensembl 116, --background-fasta would become required; with fewer, the default stays and the
warning gives the measured figure.  It changed none, so the default stays and the NOTE says so.
"""
import json

from isoform_dominance import cli

from test_fasta_competitors import CFG, SEQS, _fasta, COMP, COMPETITOR


def _run(tmp_path, capsys, *extra):
    cfg, sq = tmp_path / "cfg.json", tmp_path / "seqs.json"
    cfg.write_text(json.dumps(CFG))
    sq.write_text(json.dumps(SEQS))
    rc = cli.main(["identifiability", "--config", str(cfg), "--sequences", str(sq),
                   "--no-gene-background", *extra])
    return rc, capsys.readouterr().err


def test_without_a_fasta_the_note_gives_the_measured_figure(tmp_path, capsys):
    rc, err = _run(tmp_path, capsys)
    assert rc == 0
    note = " ".join(err.split())
    assert "not the whole index" in note
    assert "100 two-class configs at Ensembl 116" in note
    assert "changed no exit status at --min-log2fc 0.5 and one verdict" in note
    assert "Pass --background-fasta" in note


def test_with_a_fasta_there_is_no_such_note(tmp_path, capsys):
    fa = _fasta(tmp_path / "bg.fa", {COMP: COMPETITOR})
    rc, err = _run(tmp_path, capsys, "--background-fasta", fa)
    assert "not the whole index" not in " ".join(err.split())
