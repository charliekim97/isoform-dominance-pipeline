"""``scripts/parity``: the comparison of the GTF path with REST, on the GENCODE 50 extract.

The scripts take every path as an argument; on the whole GENCODE 44, 48 and 50 files they
are run by hand (``scripts/parity/README.md``).  Here they run on ``tests/data/gencode_mini``,
and ``record_rest.py`` against an offline fake of REST 116.
"""
import gzip
import importlib.util
import json
import pathlib
import urllib.request

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data" / "gencode_mini"
GTF, FASTA = str(DATA / "gencode.v50.mini.gtf.gz"), str(DATA / "gencode.v50.mini.transcripts.fa.gz")
REST = DATA / "rest116_mini.json.gz"


def _script(name):
    spec = importlib.util.spec_from_file_location(
        "parity_" + name, ROOT / "scripts" / "parity" / ("%s.py" % name))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def compare(monkeypatch):
    """compare.main, with what it patches in annotate put back afterwards."""
    from isoform_dominance import annotate
    for name in ("_get", "_post"):
        monkeypatch.setattr(annotate, name, getattr(annotate, name))
    monkeypatch.setattr(annotate.ensembl, "resolve_server", annotate.ensembl.resolve_server)
    return _script("compare").main


def _rest(tmp_path, change):
    doc = json.loads(gzip.decompress(REST.read_bytes()))
    change(doc)
    p = tmp_path / "rest.json"
    p.write_text(json.dumps(doc))
    return str(p)


@pytest.mark.parametrize("whole", [False, True])
def test_the_extract_holds_rest_116(compare, capsys, whole):
    rc = compare(["--gtf", GTF, "--fasta", FASTA, "--rest", str(REST)]
                 + (["--whole"] if whole else []))
    out = capsys.readouterr().out
    assert rc == 0, out
    assert "release 116: configs equal 6/6 (annotation_source aside); cDNA byte-identical " \
           "217/217" in out


def test_a_config_that_differs_is_named_and_exits_1(compare, capsys, tmp_path):
    def other_canonical(doc):
        g = doc["lookups"]["ENSG00000116678"]
        g["canonical_transcript"] = "ENST00000371059.8"
        for t in g["Transcript"]:
            t["is_canonical"] = int(t["id"] == "ENST00000371059")
    rc = compare(["--gtf", GTF, "--fasta", FASTA, "--rest", _rest(tmp_path, other_canonical),
                  "--json", str(tmp_path / "rows.json")])
    out = capsys.readouterr().out
    assert rc == 1
    assert "LEPR       config DIFFERS" in out and "configs equal 5/6" in out
    rows = {r["gene"]: r for r in json.loads((tmp_path / "rows.json").read_text())["rows"]}
    assert rows["LEPR"]["config_equal"] is False and rows["FOXO1"]["config_equal"] is True


def test_a_sequence_that_differs_is_named_and_exits_1(compare, capsys, tmp_path):
    def other_cdna(doc):
        doc["cdna_md5"]["ENST00000349533"] = "0" * 32
    rc = compare(["--gtf", GTF, "--fasta", FASTA, "--rest", _rest(tmp_path, other_cdna)])
    out = capsys.readouterr().out
    assert rc == 1 and "1 DIFFER: ENST00000349533" in out and "216/217" in out


def test_record_rest_writes_what_compare_reads(compare, capsys, tmp_path, monkeypatch):
    from test_gtf_live_workflow import LiveRest
    monkeypatch.setattr(urllib.request, "urlopen", LiveRest())
    genes = tmp_path / "genes.txt"
    genes.write_text("LEPR\nCD99   # a pseudoautosomal pair\n\nFOXO1\n")
    out = tmp_path / "rest116.json.gz"
    assert _script("record_rest").main(["--release", "116", "--genes", str(genes),
                                        "--out", str(out)]) == 0
    doc = json.loads(gzip.decompress(out.read_bytes()))
    assert doc["release"] == 116
    assert doc["xrefs"]["CD99"] == ["ENSG00000002586", "ENSG00000292348"]
    assert doc["lookup_symbol"]["CD99"] == "ENSG00000292348"     # REST's choice: chrY
    # every transcript of every gene, CD99's two included
    assert len(doc["cdna_md5"]) == 18 + 66 + 66 + 10
    capsys.readouterr()
    rc = compare(["--gtf", GTF, "--fasta", FASTA, "--rest", str(out)])
    assert rc == 0 and "configs equal 3/3" in capsys.readouterr().out
