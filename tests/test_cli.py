"""CLI: argument parsing, exit codes, and an offline end-to-end run."""
import json, os, random
import pytest
from isoform_dominance import cli, _selftest, annotate, identifiability


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


# transcript-scale sequences, so the default read/fragment model in
# `identifiability` is applied to something it can sensibly describe
_SHARED = _seq(1500, 21)
_ALT_A = _seq(800, 22)
_ALT_B = _seq(1000, 23)


def test_version_exits_zero(capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["--version"])
    assert e.value.code == 0
    assert "isoform-dominance" in capsys.readouterr().out


def test_no_subcommand_is_error():
    with pytest.raises(SystemExit) as e:
        cli.main([])
    assert e.value.code != 0


def _write_case(tmp_path, groups, primary, seqs):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"groups": groups, "primary_comparison": primary}))
    sq = tmp_path / "seqs.json"
    sq.write_text(json.dumps(seqs))
    return str(cfg), str(sq)


def test_identifiability_cli_ok(tmp_path):
    cfg, seqs = _write_case(
        tmp_path, {"A": ["A1"], "B": ["B1"]}, ["A", "B"],
        {"A1": _SHARED + _ALT_A, "B1": _SHARED + _ALT_B})
    rc = cli.main(["identifiability", "--config", cfg, "--sequences", seqs])
    assert rc == cli.EXIT_OK


def test_not_identifiable_without_min_log2fc_exits_zero(tmp_path, capsys):
    """The structural verdict is a flag in the report, not the exit status.

    It moves with the annotation release -- NTRK3 is `not_identifiable` on live Ensembl
    and `identifiable` on GENCODE v44 -- so it cannot gate a program's exit code.
    """
    dup = _SHARED + _ALT_A
    cfg, seqs = _write_case(
        tmp_path, {"one": ["D1"], "two": ["D2"]}, ["one", "two"],
        {"D1": dup, "D2": dup})
    rc = cli.main(["identify", "--config", cfg, "--sequences", seqs])
    out = capsys.readouterr()
    assert "VERDICT: not_identifiable" in out.out
    assert "not the exit status" in out.err
    assert rc == cli.EXIT_OK


def test_not_identifiable_with_an_unmeetable_min_log2fc_exits_nonzero(tmp_path):
    dup = _SHARED + _ALT_A
    cfg, seqs = _write_case(
        tmp_path, {"one": ["D1"], "two": ["D2"]}, ["one", "two"],
        {"D1": dup, "D2": dup})
    rc = cli.main(["identify", "--config", cfg, "--sequences", seqs,
                   "--min-log2fc", "1.0"])
    assert rc == cli.EXIT_EFFECT_NOT_RESOLVED
    assert rc != cli.EXIT_OK


def test_the_exit_status_follows_the_requested_effect_size(tmp_path):
    """Same gene, same design: only --min-log2fc moves the exit status.

    The contrast resolves |log2FC| 0.31 here, so 1.0 is met and 0.1 is not; without
    the flag the status is 0 either way.
    """
    cfg, seqs = _write_case(
        tmp_path, {"A": ["A1"], "B": ["B1"]}, ["A", "B"],
        {"A1": _SHARED + _ALT_A, "B1": _SHARED + _ALT_B})
    base = ["identify", "--config", cfg, "--sequences", seqs]
    assert cli.main(base) == cli.EXIT_OK
    assert cli.main(base + ["--min-log2fc", "1.0"]) == cli.EXIT_OK
    assert cli.main(base + ["--min-log2fc", "0.1"]) == cli.EXIT_EFFECT_NOT_RESOLVED


def test_estimable_but_starved_exits_zero_unless_an_effect_size_is_asked_for(tmp_path):
    """Nested classes: no unique k-mer, still estimable, but short of evidence.

    v2.1 exited 2 here on the strength of the zero unique-k-mer count alone, and
    v2.2 exited 3 on the verdict. The status now answers only the question asked.
    """
    cfg, seqs = _write_case(
        tmp_path, {"sub": ["S1"], "sup": ["P1"]}, ["sub", "sup"],
        {"S1": _SHARED, "P1": _SHARED + _ALT_B})
    base = ["identify", "--config", cfg, "--sequences", seqs, "--tpm", "0.01"]
    assert cli.main(base) == cli.EXIT_OK
    assert cli.main(base + ["--min-log2fc", "1.0"]) == cli.EXIT_EFFECT_NOT_RESOLVED


def test_a_transcript_shorter_than_the_window_still_exits_2(tmp_path, capsys):
    """The one case left on exit 2: a precondition, not a judgement.

    A transcript shorter than ``window`` has no windows and an all-zero column, so the
    gene total itself is outside the row space and nothing below it is well posed.
    The control drops only the short transcript, so the assertion depends on it.
    """
    seqs = {"A1": _SHARED + _ALT_A, "A2": _seq(20, 5), "B1": _SHARED + _ALT_B}
    cfg, sq = _write_case(tmp_path, {"A": ["A1", "A2"], "B": ["B1"]}, ["A", "B"], seqs)
    assert cli.main(["identify", "--config", cfg, "--sequences", sq]) == \
        cli.EXIT_NOT_IDENTIFIABLE
    assert "gene total" in capsys.readouterr().err

    assert cli.main(["identify", "--config", cfg, "--sequences", sq, "--json"]) == \
        cli.EXIT_NOT_IDENTIFIABLE
    payload = json.loads(capsys.readouterr().out)
    assert payload["gene_total"]["estimable"] is False
    assert payload["gene_total"]["transcripts_without_windows"] == ["A2"]

    ctrl = tmp_path / "ctrl"
    ctrl.mkdir()
    cfg, sq = _write_case(ctrl, {"A": ["A1"], "B": ["B1"]}, ["A", "B"],
                          {t: s for t, s in seqs.items() if t != "A2"})
    assert cli.main(["identify", "--config", cfg, "--sequences", sq]) == cli.EXIT_OK


def test_an_inestimable_contrast_is_not_a_precondition_failure(tmp_path, capsys):
    """Identical classes are `not_identifiable`, but the gene total is estimable."""
    dup = _SHARED + _ALT_A
    cfg, seqs = _write_case(
        tmp_path, {"one": ["D1"], "two": ["D2"]}, ["one", "two"],
        {"D1": dup, "D2": dup})
    cli.main(["identify", "--config", cfg, "--sequences", seqs, "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["verdict"] == "not_identifiable"
    assert payload["gene_total"]["estimable"] is True
    assert payload["effect_resolvable"] is None


def test_identifiability_cli_json(tmp_path, capsys):
    cfg, seqs = _write_case(
        tmp_path, {"A": ["A1"], "B": ["B1"]}, ["A", "B"],
        {"A1": _SHARED + _ALT_A, "B1": _SHARED + _ALT_B})
    rc = cli.main(["identifiability", "--config", cfg, "--sequences", seqs, "--json"])
    assert rc == cli.EXIT_OK
    payload = json.loads(capsys.readouterr().out)
    assert payload["verdict"] == "identifiable"
    assert payload["groups"]["A"]["expected_informative_reads"] > 0
    assert payload["contrast"]["estimable"] is True


def test_identifiability_cli_background_fasta_changes_the_verdict(tmp_path, capsys):
    fa = tmp_path / "bg.fa"
    fa.write_text(">DECOY\n%s\n" % (_SHARED + _ALT_A))
    cfg, seqs = _write_case(
        tmp_path, {"A": ["A1"], "B": ["B1"]}, ["A", "B"],
        {"A1": _SHARED + _ALT_A, "B1": _SHARED + _ALT_B})
    base = ["identify", "--config", cfg, "--sequences", seqs, "--json"]
    cli.main(base)
    assert json.loads(capsys.readouterr().out)["verdict"] == "identifiable"
    cli.main(base + ["--background-fasta", str(fa)])
    assert json.loads(capsys.readouterr().out)["verdict"] != "identifiable"


def test_annotate_cli_offline(tmp_path, monkeypatch):
    out = tmp_path / "auto.json"

    def fake_run(gene, outpath, species="homo_sapiens", **retry):
        cfg = {"gene": gene, "groups": {"iso_100aa": ["T1"], "iso_50aa": ["T2"]},
               "primary_comparison": ["iso_50aa", "iso_100aa"]}
        json.dump(cfg, open(outpath, "w"))
        return cfg

    monkeypatch.setattr(annotate, "run", fake_run)
    rc = cli.main(["annotate", "--gene", "FAKE", "--out", str(out)])
    assert (rc or 0) == 0
    assert json.load(open(out))["gene"] == "FAKE"


@pytest.mark.parametrize("sub", ["annotate", "identifiability"])
def test_a_read_timeout_is_a_reported_network_failure_not_a_traceback(
        sub, tmp_path, monkeypatch, capsys):
    """A read timeout raises TimeoutError, which is not a URLError, so it slipped past
    `_net_fail` and reached the user as a traceback."""
    def timed_out(*a, **k):
        raise TimeoutError("The read operation timed out")

    monkeypatch.setattr(annotate, "run", timed_out)
    monkeypatch.setattr(identifiability, "analyze", timed_out)
    cfg, _ = _write_case(tmp_path, {"A": ["A1"], "B": ["B1"]}, ["A", "B"], {})
    argv = (["annotate", "--gene", "X", "--out", str(tmp_path / "o.json")]
            if sub == "annotate" else ["identifiability", "--config", cfg])
    assert cli.main(argv) == 1
    err = capsys.readouterr().err
    assert "Ensembl request failed" in err
    assert "timed out" in err


def test_extract_stats_qc_cli_end_to_end(tmp_path):
    base = str(tmp_path)
    info = _selftest.generate(base)
    cfg_path = os.path.join(base, "config.json")
    json.dump(_selftest.CONFIG, open(cfg_path, "w"))

    perdonor = {}
    for cohort, p in info.items():
        out = os.path.join(base, "pd_%s.csv" % cohort)
        rc = cli.main(["extract", "--config", cfg_path, "--quantdir", p["quantdir"],
                       "--samplemap", p["samplemap"], "--cohort", cohort, "--out", out])
        assert (rc or 0) == 0
        perdonor[cohort] = out

    rc = cli.main(["stats", "--config", cfg_path, "--condition", "control"]
                  + sum([["--perdonor", "%s=%s" % (c, perdonor[c])] for c in perdonor], [])
                  + ["--out", os.path.join(base, "res")])
    assert (rc or 0) == 0
    assert os.path.exists(os.path.join(base, "res.png"))

    rc = cli.main(["qc", "--config", cfg_path]
                  + sum([["--markers", "%s=%s" % (c, info[c]["markers"])] for c in info], [])
                  + sum([["--target", "%s=%s" % (c, perdonor[c])] for c in perdonor], [])
                  + ["--out", os.path.join(base, "qc")])
    assert (rc or 0) == 0
    assert os.path.exists(os.path.join(base, "qc.png"))
