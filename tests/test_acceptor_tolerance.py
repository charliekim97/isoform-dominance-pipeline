"""``annotate --acceptor-tolerance N`` (issue #4): transcripts whose terminal-exon acceptors
lie within N bp of a cluster's first acceptor are one class.

The rule: acceptors in ascending order; a cluster takes every acceptor within N bp of its
first (lowest) one and no further, so N bounds a cluster's span and nothing chains; the
cluster's ``terminal_acceptor`` is that first coordinate.  N = 0, the default, is the exact
coordinate, and its config is 2.6.0's byte for byte (``test_pinned_outputs.py``).  REST and
a GTF go through the same function, so at every N they propose alike.
"""
import gzip
import json
import pathlib

import pytest

from isoform_dominance import annotate, cli

DATA = pathlib.Path(__file__).parent / "data"
GTF = DATA / "gencode_mini" / "gencode.v50.mini.gtf.gz"
GOLDEN = DATA / "v260_outputs"
REST = json.loads(gzip.decompress((DATA / "gencode_mini" / "rest116_mini.json.gz")
                                  .read_bytes()))
BY_NAME = {}
for _g in REST["lookups"].values():
    BY_NAME.setdefault(_g["display_name"], []).append(_g)
GENES = ["LEPR", "FOXO1", "STK11", "AXIN1", "GSK3B", "CD99"]


@pytest.fixture
def rest(monkeypatch):
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


def _info(accs):
    """One transcript per acceptor, protein length 100 + its index."""
    return {"transcripts": [{"id": "T%d" % i, "protein_aa": 100 + i, "terminal_acceptor": a,
                             "is_canonical": i == 0} for i, a in enumerate(accs)]}


def _clusters(accs, n):
    return sorted(sorted(int(t[1:]) for t in c["ids"])
                  for c in annotate.cluster_by_terminal_exon(_info(accs), tolerance=n))


# --------------------------------------------------------------------------- #
# the rule
# --------------------------------------------------------------------------- #
def test_a_cluster_takes_what_lies_within_n_of_its_first_acceptor_and_no_further():
    accs = [100, 103, 106]
    assert _clusters(accs, 0) == [[0], [1], [2]]
    assert _clusters(accs, 2) == [[0], [1], [2]]
    assert _clusters(accs, 3) == [[0, 1], [2]]         # 106 is 6 from 100: no chaining
    assert _clusters(accs, 6) == [[0, 1, 2]]


def test_the_order_the_transcripts_come_in_does_not_matter():
    accs = [106, 100, 103, 230, 228]
    for n in (0, 2, 3, 6):
        got = annotate.cluster_by_terminal_exon(_info(accs), tolerance=n)
        back = annotate.cluster_by_terminal_exon(
            {"transcripts": list(reversed(_info(accs)["transcripts"]))}, tolerance=n)
        assert got == back


def test_a_merged_cluster_is_at_its_first_acceptor_and_lists_its_members():
    got = annotate.cluster_by_terminal_exon(_info([500, 497, 499, 497]), tolerance=3)
    assert len(got) == 1
    c = got[0]
    assert c["acceptor"] == 497
    assert c["ids"] == ["T0", "T1", "T2", "T3"]
    assert c["members"] == [{"acceptor": 497, "offset": 0, "ids": ["T1", "T3"]},
                            {"acceptor": 499, "offset": 2, "ids": ["T2"]},
                            {"acceptor": 500, "offset": 3, "ids": ["T0"]}]
    assert c["span"] == 3
    # a cluster of one acceptor is what 2.6.0 wrote, at any tolerance
    alone = annotate.cluster_by_terminal_exon(_info([10, 900]), tolerance=3)
    assert all(set(c) == {"acceptor", "rep_aa", "n", "canonical", "ids"} for c in alone)


def test_a_negative_tolerance_is_refused():
    with pytest.raises(ValueError, match="tolerance"):
        annotate.cluster_by_terminal_exon(_info([1, 2]), tolerance=-1)


# --------------------------------------------------------------------------- #
# REST and a GTF: one function
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("n", [3, 6, 10, 20])
@pytest.mark.parametrize("gene", GENES)
def test_the_two_sources_propose_alike_at_every_tolerance(rest, gene, n):
    by_file = annotate.build_config_from_gtf(gene, str(GTF), acceptor_tolerance=n)
    by_file.pop("annotation_source")
    assert by_file == annotate.build_config(gene, acceptor_tolerance=n)


def test_on_the_extract_a_tolerance_of_six_merges_acceptors_of_gsk3b_and_cd99(rest):
    """The parity above is not vacuous: GSK3B has two acceptors 4 bp apart, CD99 two 5 bp
    apart, and the other four genes none within 20 bp."""
    for gene in GENES:
        cfg = annotate.build_config_from_gtf(gene, str(GTF), acceptor_tolerance=6)
        merged = [c for c in cfg["_clusters"] if "merged_acceptors" in c]
        assert bool(merged) == (gene in ("GSK3B", "CD99")), gene
        for c in merged:
            assert c["acceptor_span"] <= 6
            assert c["terminal_acceptor"] == c["merged_acceptors"][0]["terminal_acceptor"]


# --------------------------------------------------------------------------- #
# the command
# --------------------------------------------------------------------------- #
def test_the_option_reaches_both_modes(rest, tmp_path, capsys):
    for mode in ([], ["--gtf", str(GTF)]):
        out = tmp_path / "cfg.json"
        assert cli.main(["annotate", "--gene", "CD99", "--acceptor-tolerance", "6",
                         "--out", str(out)] + mode) == 0
        cfg = json.loads(out.read_text())
        assert cfg["_proposal"]["acceptor_tolerance"] == 6
        assert any("merged_acceptors" in c for c in cfg["_clusters"])
        assert "6 bp" in cfg["_proposed"]
    capsys.readouterr()


@pytest.mark.parametrize("mode", ["rest", "gtf"])
def test_tolerance_zero_given_outright_is_the_2_6_0_config(rest, tmp_path, mode):
    out = tmp_path / "cfg.json"
    extra = ["--gtf", str(GTF)] if mode == "gtf" else []
    assert cli.main(["annotate", "--gene", "GSK3B", "--acceptor-tolerance", "0",
                     "--out", str(out)] + extra) == 0
    golden = (GOLDEN / ("annotate_%s_GSK3B.json" % mode)).read_text()
    assert out.read_text() == golden


@pytest.mark.parametrize("bad", ["-1", "x", "1.5"])
def test_the_command_refuses_a_tolerance_that_is_not_a_count(tmp_path, capsys, bad):
    with pytest.raises(SystemExit) as e:
        cli.main(["annotate", "--gene", "LEPR", "--gtf", str(GTF), "--acceptor-tolerance",
                  bad, "--out", str(tmp_path / "x.json")])
    assert e.value.code == 2
    assert "--acceptor-tolerance" in capsys.readouterr().err
