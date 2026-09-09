"""CLI: argument parsing, exit codes, and an offline end-to-end run."""
import json, os, random
import pytest
from isoform_dominance import cli, _selftest, annotate


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


def test_identifiability_cli_exit2_when_classes_are_indistinguishable(tmp_path):
    # identical sequence: neither class total nor their contrast is estimable,
    # and no amount of sequencing changes that
    dup = _SHARED + _ALT_A
    cfg, seqs = _write_case(
        tmp_path, {"one": ["D1"], "two": ["D2"]}, ["one", "two"],
        {"D1": dup, "D2": dup})
    rc = cli.main(["identify", "--config", cfg, "--sequences", seqs])
    assert rc == cli.EXIT_NOT_IDENTIFIABLE


def test_identifiability_cli_exit3_when_estimable_but_starved(tmp_path):
    """Nested classes: no unique k-mer, still estimable, but short of evidence.

    v2.1 exited 2 here on the strength of the zero unique-k-mer count alone.  The
    contrast is in fact recoverable, so the honest answer is the middle one.
    """
    cfg, seqs = _write_case(
        tmp_path, {"sub": ["S1"], "sup": ["P1"]}, ["sub", "sup"],
        {"S1": _SHARED, "P1": _SHARED + _ALT_B})
    rc = cli.main(["identify", "--config", cfg, "--sequences", seqs,
                   "--tpm", "0.01"])
    assert rc == cli.EXIT_WEAK


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


def test_identifiability_cli_background_fasta_changes_the_verdict(tmp_path):
    fa = tmp_path / "bg.fa"
    fa.write_text(">DECOY\n%s\n" % (_SHARED + _ALT_A))
    cfg, seqs = _write_case(
        tmp_path, {"A": ["A1"], "B": ["B1"]}, ["A", "B"],
        {"A1": _SHARED + _ALT_A, "B1": _SHARED + _ALT_B})
    assert cli.main(["identify", "--config", cfg, "--sequences", seqs]) == cli.EXIT_OK
    assert cli.main(["identify", "--config", cfg, "--sequences", seqs,
                     "--background-fasta", str(fa)]) != cli.EXIT_OK


def test_annotate_cli_offline(tmp_path, monkeypatch):
    out = tmp_path / "auto.json"

    def fake_run(gene, outpath, species="homo_sapiens"):
        cfg = {"gene": gene, "groups": {"iso_100aa": ["T1"], "iso_50aa": ["T2"]},
               "primary_comparison": ["iso_50aa", "iso_100aa"]}
        json.dump(cfg, open(outpath, "w"))
        return cfg

    monkeypatch.setattr(annotate, "run", fake_run)
    rc = cli.main(["annotate", "--gene", "FAKE", "--out", str(out)])
    assert (rc or 0) == 0
    assert json.load(open(out))["gene"] == "FAKE"


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
