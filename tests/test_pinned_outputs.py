"""What 2.7.0 must not move: the self-test's reference P values, `annotate`'s config at the
default acceptor tolerance, and `extract`'s output.

The files compared against are in ``tests/data/v260_outputs`` (its README says how they
were written).  JSON is written in text mode, so on Windows its newlines are ``\\r\\n`` --
as 2.6.0's were there -- and the golden bytes are compared in the platform's newline.
"""
import gzip
import json
import os
import pathlib

import pytest

from isoform_dominance import _selftest, annotate, cli, extract, stats

DATA = pathlib.Path(__file__).parent / "data"
GOLDEN = DATA / "v260_outputs"
MINI = DATA / "gencode_mini"
GTF = MINI / "gencode.v50.mini.gtf.gz"
REST = json.loads(gzip.decompress((MINI / "rest116_mini.json.gz").read_bytes()))
BY_NAME = {}
for _g in REST["lookups"].values():
    BY_NAME.setdefault(_g["display_name"], []).append(_g)
GENES = ["LEPR", "FOXO1", "STK11", "AXIN1", "GSK3B", "CD99"]


def _golden_text(name):
    """A golden JSON file's bytes in this platform's newline, as text mode writes them."""
    return (GOLDEN / name).read_bytes().replace(b"\n", os.linesep.encode())


@pytest.fixture
def rest(monkeypatch):
    """annotate's REST calls answered from the recorded lookups (as in test_gtf_parity)."""
    def get(path, **kw):
        if path.startswith("/info/data"):
            return {"releases": [REST["release"]]}
        sym = path.split("?")[0].split("/")[4] if "/symbol/" in path else None
        if path.startswith("/lookup/symbol/"):
            gid = REST["lookup_symbol"].get(sym) or BY_NAME[sym][0]["id"]
            return REST["lookups"][gid]
        if path.startswith("/xrefs/symbol/"):
            return [{"type": "gene", "id": g["id"]} for g in BY_NAME.get(sym, [])]
        if path.startswith("/lookup/id/"):
            return REST["lookups"][path.split("/")[3].split("?")[0]]
        raise AssertionError(path)
    monkeypatch.setattr(annotate, "_get", get)
    monkeypatch.setattr(annotate, "_post", lambda path, body, **kw:
                        {i: REST["lookups"][i] for i in body["ids"] if i in REST["lookups"]})
    monkeypatch.setattr(annotate.ensembl, "resolve_server", lambda release=None, **kw: None)


# --------------------------------------------------------------------------- #
# the self-test's reference result
# --------------------------------------------------------------------------- #
def _selftest_stats(tmp_path):
    info = _selftest.generate(str(tmp_path))
    perdonor = {}
    for cohort, p in info.items():
        out = str(tmp_path / ("pd_%s.csv" % cohort))
        extract.run(_selftest.CONFIG, p["quantdir"], p["samplemap"], cohort, out)
        perdonor[cohort] = out
    return stats.run(_selftest.CONFIG, "control", perdonor, str(tmp_path / "r"))


def test_the_self_tests_pooled_and_stouffer_p_do_not_move(tmp_path):
    """Both cohorts are unanimous with no zero difference, so neither the Stouffer sign
    rule nor its weights can move the combination; 2.6.0's values, to the last digit."""
    res = _selftest_stats(tmp_path)
    cn, cgt, cp, _ = res["combined"]
    assert (cn, cgt) == (11, 11)
    assert cp == 2 ** -10                                   # 9.766e-4
    assert "%.4g" % cp == "0.0009766"
    stouffer = res["combination"]["stouffer"]
    assert stouffer["p"] == pytest.approx(0.004418947896206985, rel=1e-12)
    assert "%.4g" % stouffer["p"] == "0.004419"
    assert stouffer["z"] == pytest.approx(2.8465954536696483, rel=1e-12)


# --------------------------------------------------------------------------- #
# annotate at the default tolerance: 2.6.0's config, byte for byte
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("gene", GENES)
def test_annotate_from_rest_writes_the_2_6_0_config(rest, tmp_path, gene):
    out = tmp_path / "cfg.json"
    assert cli.main(["annotate", "--gene", gene, "--out", str(out)]) == 0
    assert out.read_bytes() == _golden_text("annotate_rest_%s.json" % gene)


@pytest.mark.parametrize("gene", GENES)
def test_annotate_from_a_gtf_writes_the_2_6_0_config(tmp_path, gene):
    out = tmp_path / "cfg.json"
    assert cli.main(["annotate", "--gene", gene, "--gtf", str(GTF), "--out", str(out)]) == 0
    assert out.read_bytes() == _golden_text("annotate_gtf_%s.json" % gene)


# --------------------------------------------------------------------------- #
# extract: untouched
# --------------------------------------------------------------------------- #
def test_extract_writes_what_2_6_0_wrote(tmp_path):
    info = _selftest.generate(str(tmp_path / "work"))
    for cohort, p in info.items():
        out = tmp_path / ("perdonor_%s.csv" % cohort)
        assert cli.main(["extract", "--config", _config_file(tmp_path),
                         "--quantdir", p["quantdir"], "--samplemap", p["samplemap"],
                         "--cohort", cohort, "--out", str(out)]) == 0
        # csv.writer on a newline="" file: the same bytes on every platform
        assert out.read_bytes() == (GOLDEN / ("extract_selftest_%s.csv" % cohort)).read_bytes()
        assert pathlib.Path(str(out) + ".index.json").read_bytes() == _golden_text(
            "extract_selftest_%s.csv.index.json" % cohort)


def _config_file(tmp_path):
    path = tmp_path / "config.json"
    if not path.exists():
        path.write_text(json.dumps(_selftest.CONFIG))
    return str(path)


def test_the_golden_files_are_the_ones_the_readme_lists():
    names = sorted(p.name for p in GOLDEN.iterdir() if p.name != "README.md")
    assert names == sorted(["annotate_rest_%s.json" % g for g in GENES]
                           + ["annotate_gtf_%s.json" % g for g in GENES]
                           + ["extract_selftest_%s.csv%s" % (c, x) for c in _selftest.DATA
                              for x in ("", ".index.json")]
                           + ["identifiability_inputs_GENEX.json"])
