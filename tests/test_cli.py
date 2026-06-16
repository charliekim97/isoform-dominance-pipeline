"""CLI: argument parsing, exit codes, and an offline end-to-end run."""
import json, os
import pytest
from isoform_dominance import cli, _selftest, annotate


def test_version_exits_zero(capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["--version"])
    assert e.value.code == 0
    assert "isoform-dominance" in capsys.readouterr().out


def test_no_subcommand_is_error():
    with pytest.raises(SystemExit) as e:
        cli.main([])
    assert e.value.code != 0


def test_identifiability_cli_ok(tmp_path):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"groups": {"A": ["A1"], "B": ["B1"]},
                               "primary_comparison": ["A", "B"]}))
    seqs = tmp_path / "seqs.json"
    seqs.write_text(json.dumps({"A1": "ACGT" * 20 + "GGGGGGGGCATCAT",
                                "B1": "ACGT" * 20 + "TTTTTTTTAGAGAG"}))
    rc = cli.main(["identifiability", "--config", str(cfg),
                   "--sequences", str(seqs), "--k", "8"])
    assert rc == 0


def test_identifiability_cli_not_distinguishable_exit2(tmp_path):
    sub = "ACGT" * 20 + "GGGGGGGG"
    sup = sub + "TTTTTTTTTT"
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"groups": {"sub": ["S1"], "super": ["P1"]},
                               "primary_comparison": ["sub", "super"]}))
    seqs = tmp_path / "seqs.json"
    seqs.write_text(json.dumps({"S1": sub, "P1": sup}))
    rc = cli.main(["identify", "--config", str(cfg),
                   "--sequences", str(seqs), "--k", "8"])
    assert rc == 2


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
