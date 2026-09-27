"""`--save-inputs` / `--inputs`: a verdict that outlives Ensembl's REST archive.

Ensembl keeps REST archives for recent releases only, and one can be down for a day:
release 111's timed out on every lookup on 2026-09-24.  A saved-inputs file holds the
sequence a run used, with the release it came from, and repeats the run with no request.
"""
import json
import time
import urllib.request

import pytest

from isoform_dominance import cli, identifiability, io
from test_ensembl_release import (BACKGROUND_AT, CONFIG, CURRENT, GROUP_A, GROUP_B, OLD,
                                  ONLY_NOW, SEQS, Ensembl)


@pytest.fixture
def ens(monkeypatch):
    monkeypatch.setattr(time, "sleep", lambda s: None)
    fake = Ensembl()
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    return fake


@pytest.fixture
def no_network(monkeypatch):
    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)


def _captured(sequences, background=None, k=31, window=None, canonical=True):
    """What ``analyze(..., inputs_out=...)`` hands back, for a direct `io.save_inputs`."""
    return {"sequences": sequences, "background_sequences": background or {},
            "k": k, "window": window if window is not None else k, "canonical": canonical}


def _cfg(tmp_path, **kw):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(dict(CONFIG, **kw)))
    return str(p)


def _same_answer(a, b):
    for key in ("verdict", "reasons", "effect_resolvable", "n_compatibility_classes"):
        assert a[key] == b[key], key
    for g in ("A", "B"):
        for key in ("n_unique_kmers", "conditioning_factor", "min_resolvable_log2fc"):
            assert a["groups"][g][key] == b["groups"][g][key], (g, key)
    assert a["contrast"]["min_resolvable_log2fc"] == b["contrast"]["min_resolvable_log2fc"]


def test_analyze_hands_back_exactly_what_it_used(ens):
    got = {}
    identifiability.analyze(CONFIG, ensembl_release=OLD, inputs_out=got)
    assert got["fetched_release"] == OLD
    assert got["sequences"] == {t: SEQS[OLD][t] for t in GROUP_A + GROUP_B}
    assert got["background_sequences"] == {t: SEQS[OLD][t] for t in BACKGROUND_AT[OLD]}


def test_a_saved_run_repeats_offline_with_the_same_answer(ens, tmp_path, capsys, monkeypatch):
    cfg, saved = _cfg(tmp_path, ensembl_release=OLD), str(tmp_path / "run.inputs.json")
    assert cli.main(["identify", "--config", cfg, "--ensembl-release", str(OLD),
                     "--save-inputs", saved, "--json"]) == 0
    online = json.loads(capsys.readouterr().out)
    doc = json.loads(open(saved).read())
    assert doc["format"] == io.INPUTS_FORMAT and doc["ensembl_release"] == OLD
    assert doc["config_ensembl_release"] == OLD and doc["gene"] == "FAKE"

    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    assert cli.main(["identify", "--config", cfg, "--inputs", saved, "--json"]) == 0
    offline = json.loads(capsys.readouterr().out)
    _same_answer(online, offline)
    assert offline["annotation"] == {"ensembl_release": OLD, "fetched_release": None,
                                     "inputs_release": OLD}


def test_the_header_names_the_saved_inputs(no_network, tmp_path, capsys):
    saved = tmp_path / "s.json"
    io.save_inputs(str(saved), _captured({t: SEQS[OLD][t] for t in GROUP_A + GROUP_B}),
                   dict(CONFIG), OLD, "test")
    cli.main(["identify", "--config", _cfg(tmp_path, ensembl_release=CURRENT),
              "--inputs", str(saved)])
    out = capsys.readouterr()
    assert "sequence from saved inputs (release %d" % OLD in out.out
    assert ("config was annotated against Ensembl release %d, and the saved inputs are "
            "from release %d" % (CURRENT, OLD)) in out.err


def test_inputs_missing_a_configured_transcript_are_refused_not_topped_up(no_network,
                                                                          tmp_path, capsys):
    saved = tmp_path / "s.json"
    io.save_inputs(str(saved), _captured({t: SEQS[CURRENT][t] for t in GROUP_A + GROUP_B}),
                   dict(CONFIG), CURRENT, "test")
    cfg = _cfg(tmp_path, groups={"A": GROUP_A, "B": GROUP_B + [ONLY_NOW]})
    assert cli.main(["identify", "--config", cfg, "--inputs", str(saved)]) == 1
    assert "saved inputs have no sequence for %s" % ONLY_NOW in capsys.readouterr().err


@pytest.mark.parametrize("extra", [["--ensembl-release", "110"], ["--sequences", "x.json"]])
def test_inputs_replace_the_other_sequence_sources(no_network, tmp_path, capsys, extra):
    saved = tmp_path / "s.json"
    io.save_inputs(str(saved), _captured({}), dict(CONFIG), OLD, "test")
    assert cli.main(["identify", "--config", _cfg(tmp_path), "--inputs", str(saved)]
                    + extra) == 1
    assert "--inputs replaces" in capsys.readouterr().err


def test_a_file_that_is_not_saved_inputs_is_refused(no_network, tmp_path, capsys):
    bogus = tmp_path / "seqs.json"
    bogus.write_text(json.dumps({t: SEQS[OLD][t] for t in GROUP_A}))   # a --sequences file
    assert cli.main(["identify", "--config", _cfg(tmp_path), "--inputs", str(bogus)]) == 1
    assert "not a saved-inputs file" in capsys.readouterr().err


def test_inputs_for_another_gene_are_refused(no_network, tmp_path, capsys):
    saved = tmp_path / "s.json"
    io.save_inputs(str(saved), _captured({t: SEQS[OLD][t] for t in GROUP_A + GROUP_B}),
                   dict(CONFIG, gene="OTHER"), OLD, "test")
    assert cli.main(["identify", "--config", _cfg(tmp_path), "--inputs", str(saved)]) == 1
    assert "for OTHER and the config for FAKE" in capsys.readouterr().err


def _save(tmp_path, release=OLD, config=None, fasta=None):
    saved = tmp_path / "s.json"
    io.save_inputs(str(saved), _captured({t: SEQS[OLD][t] for t in GROUP_A + GROUP_B}),
                   dict(config or CONFIG), release, "test", background_fasta=fasta)
    return str(saved)


def test_no_release_is_named_none(no_network, tmp_path, capsys):
    # a config without a release, and inputs without one: say so, never "release None"
    cfg = {k: v for k, v in CONFIG.items() if k != "ensembl_release"}
    p = tmp_path / "c.json"
    p.write_text(json.dumps(cfg))
    cli.main(["identify", "--config", str(p), "--inputs", _save(tmp_path)])
    out = capsys.readouterr()
    assert "None" not in out.out + out.err
    cli.main(["identify", "--config", _cfg(tmp_path, ensembl_release=CURRENT),
              "--inputs", _save(tmp_path, release=None)])
    out = capsys.readouterr()
    assert "None" not in out.out + out.err
    assert "saved inputs (release not recorded" in out.out
    assert "record no release: their sequence was supplied offline" in out.err


def _fasta(tmp_path, name, seq):
    p = tmp_path / name
    p.write_text(">ENSTOTHER0001 other gene\n%s\n" % seq)
    return str(p)


def test_the_background_fasta_is_recorded_and_its_absence_or_change_is_named(
        no_network, tmp_path, capsys):
    # the FASTA shares a stretch of group A, so it changes the verdict: a rerun that
    # silently dropped or swapped it would not be the run that was saved
    shared = SEQS[OLD][GROUP_A[0]][100:400]
    fa = _fasta(tmp_path, "bg.fa", shared)
    cfg = _cfg(tmp_path, ensembl_release=OLD)
    saved = _save(tmp_path, fasta=fa)
    doc = json.loads(open(saved).read())
    assert doc["background_fasta"]["path"] == fa
    assert doc["background_fasta"]["bytes"] == len(open(fa, "rb").read())
    assert doc["background_fasta"]["sha256"] == io.file_sha256(fa)

    def run(*extra):
        code = cli.main(["identify", "--config", cfg, "--inputs", saved, "--json"]
                        + list(extra))
        out = capsys.readouterr()
        return code, json.loads(out.out), out.err

    _, same, err = run("--background-fasta", fa)
    assert "NOTE" not in err
    _, without, err = run()
    assert "also used --background-fasta %s" % fa in err
    assert same["groups"]["A"]["n_unique_kmers"] < without["groups"]["A"]["n_unique_kmers"]
    other = _fasta(tmp_path, "other.fa", shared[:-1] + ("A" if shared[-1] != "A" else "C"))
    _, _, err = run("--background-fasta", other)
    assert "--background-fasta %s is not the file the inputs were saved with" % other in err


def test_a_run_saved_without_a_fasta_says_nothing_about_one(no_network, tmp_path, capsys):
    assert json.loads(open(_save(tmp_path)).read())["background_fasta"] is None
    cli.main(["identify", "--config", _cfg(tmp_path, ensembl_release=OLD),
              "--inputs", _save(tmp_path)])
    err = capsys.readouterr().err
    assert "also used --background-fasta" not in err and "inputs were saved with" not in err


def test_an_unwritable_save_path_is_refused_before_anything_is_fetched(tmp_path, capsys,
                                                                       monkeypatch):
    def refuse(req, *a, **k):
        raise AssertionError("fetched before the save path was checked: %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    bad = tmp_path / "no-such-dir" / "run.json"
    assert cli.main(["identify", "--config", _cfg(tmp_path, ensembl_release=OLD),
                     "--ensembl-release", str(OLD), "--save-inputs", str(bad)]) == 1
    err = capsys.readouterr().err
    assert "--save-inputs %s: directory" % bad in err and len(err.strip().splitlines()) == 1


def test_a_save_path_that_is_a_directory_is_refused_before_anything_is_fetched(
        tmp_path, capsys, monkeypatch):
    # dirname() of a directory is its parent, which exists, so this used to pass the
    # check and fail on open() -- after the run, with the fetched sequence thrown away
    def refuse(req, *a, **k):
        raise AssertionError("fetched before the save path was checked: %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    already = tmp_path / "already-a-directory"
    already.mkdir()
    assert cli.main(["identify", "--config", _cfg(tmp_path, ensembl_release=OLD),
                     "--ensembl-release", str(OLD), "--save-inputs", str(already)]) == 1
    err = capsys.readouterr().err
    assert "--save-inputs %s: is a directory" % already in err
    assert len(err.strip().splitlines()) == 1


# --------------------------------------------------------------------------- #
# the config decides which saved sequence is a class and which is background
# --------------------------------------------------------------------------- #
def _saved_run(tmp_path, capsys):
    """Save a run of the whole config at release OLD; return (saved path, its doc)."""
    saved = str(tmp_path / "run.inputs.json")
    assert cli.main(["identify", "--config", _cfg(tmp_path, ensembl_release=OLD),
                     "--ensembl-release", str(OLD), "--save-inputs", saved, "--json"]) == 0
    capsys.readouterr()
    return saved


def _run(argv, capsys, expect=0):
    assert cli.main(argv + ["--json"]) == expect
    return json.loads(capsys.readouterr().out)


def test_a_transcript_the_config_no_longer_names_becomes_background(ens, tmp_path, capsys,
                                                                    monkeypatch):
    # the saved file holds its sequence; dropping it from the groups must put it where a
    # live run puts it -- in the gene background -- not throw it away
    saved = _saved_run(tmp_path, capsys)
    dropped = GROUP_A[-1]
    cfg = _cfg(tmp_path, ensembl_release=OLD, groups={"A": GROUP_A[:-1], "B": GROUP_B})
    live = _run(["identify", "--config", cfg, "--ensembl-release", str(OLD)], capsys)
    assert dropped in live["background"]["gene_transcripts"]

    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    offline = _run(["identify", "--config", cfg, "--inputs", saved], capsys)
    assert dropped in offline["background"]["gene_transcripts"]
    assert (offline["background"]["n_background_transcripts"]
            == live["background"]["n_background_transcripts"])
    _same_answer(live, offline)


def test_a_transcript_moved_from_the_background_into_a_class_is_accepted(ens, tmp_path,
                                                                         capsys, monkeypatch):
    # its sequence is in the file, under background_sequences: a rerun must use it, not
    # refuse the config for a sequence it is holding
    saved = _saved_run(tmp_path, capsys)
    promoted = BACKGROUND_AT[OLD][0]
    cfg = _cfg(tmp_path, ensembl_release=OLD,
               groups={"A": GROUP_A, "B": GROUP_B + [promoted]})
    live = _run(["identify", "--config", cfg, "--ensembl-release", str(OLD)], capsys)
    assert promoted not in live["background"]["gene_transcripts"]

    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    offline = _run(["identify", "--config", cfg, "--inputs", saved], capsys)
    assert promoted not in offline["background"]["gene_transcripts"]
    assert offline["groups"]["B"]["n_transcripts"] == len(GROUP_B) + 1
    _same_answer(live, offline)


def test_a_transcript_no_release_ever_had_is_still_refused(ens, tmp_path, capsys,
                                                           monkeypatch):
    # the merge above must not turn the "saved for another grouping" refusal vacuous
    saved = _saved_run(tmp_path, capsys)
    cfg = _cfg(tmp_path, ensembl_release=OLD,
               groups={"A": GROUP_A, "B": GROUP_B + ["ENST00000999999"]})

    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    assert cli.main(["identify", "--config", cfg, "--inputs", saved]) == 1
    assert "saved inputs have no sequence for ENST00000999999" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# the file records the analysis it was saved from
# --------------------------------------------------------------------------- #
def test_the_saved_file_records_the_analysis_parameters(ens, tmp_path, capsys):
    saved = str(tmp_path / "run.inputs.json")
    assert cli.main(["identify", "--config", _cfg(tmp_path, ensembl_release=OLD),
                     "--ensembl-release", str(OLD), "--save-inputs", saved,
                     "--k", "25", "--window", "40", "--strand-aware", "--json"]) == 0
    capsys.readouterr()
    doc = json.loads(open(saved).read())
    assert doc["analysis"] == {"k": 25, "window": 40, "canonical": False}


def test_a_rerun_at_another_window_or_k_says_so(ens, tmp_path, capsys, monkeypatch):
    # k, the window and the k-mer convention are not recorded in the config, so nothing
    # else can tell a rerun that it is not repeating the run that was saved
    saved = str(tmp_path / "run.inputs.json")
    cfg = _cfg(tmp_path, ensembl_release=OLD)
    assert cli.main(["identify", "--config", cfg, "--ensembl-release", str(OLD),
                     "--save-inputs", saved, "--k", "25", "--json"]) == 0
    capsys.readouterr()

    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    assert cli.main(["identify", "--config", cfg, "--inputs", saved, "--k", "25"]) == 0
    assert "the inputs were saved" not in capsys.readouterr().err
    assert cli.main(["identify", "--config", cfg, "--inputs", saved, "--k", "31"]) == 0
    assert ("the inputs were saved at k=25, window=25, canonical k-mers; this run is at "
            "k=31, window=31, canonical k-mers") in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# a malformed saved-inputs file is one line, not a traceback
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("drop,expect", [
    ("ensembl_release", "records no ensembl_release"),
    ("analysis", "records no analysis"),
])
def test_a_saved_inputs_file_missing_a_key_is_one_line(no_network, tmp_path, capsys,
                                                       drop, expect):
    saved = tmp_path / "s.json"
    doc = json.loads(open(_save(tmp_path)).read())
    del doc[drop]
    saved.write_text(json.dumps(doc))
    assert cli.main(["identify", "--config", _cfg(tmp_path), "--inputs", str(saved)]) == 1
    err = capsys.readouterr().err
    assert expect in err and len(err.strip().splitlines()) == 1


@pytest.mark.parametrize("bad,field", [
    ({"sequences": {"ENST00000000001": 17}}, "sequences"),
    ({"background_sequences": [1, 2, 3]}, "background_sequences"),
    ({"background_fasta": {"path": "x.fa"}}, "background_fasta"),
    ({"analysis": {"k": "thirty-one"}}, "analysis"),
    ({"ensembl_release": "110"}, "ensembl_release"),
])
def test_a_saved_inputs_file_with_a_bad_value_is_one_line(no_network, tmp_path, capsys,
                                                          bad, field):
    saved = tmp_path / "s.json"
    saved.write_text(json.dumps(dict(json.loads(open(_save(tmp_path)).read()), **bad)))
    assert cli.main(["identify", "--config", _cfg(tmp_path), "--inputs", str(saved)]) == 1
    err = capsys.readouterr().err
    assert field in err and str(saved) in err
    assert len(err.strip().splitlines()) == 1 and "Traceback" not in err
