"""What `identifiability` says when Ensembl does not have what the config names.

A ``gene_id`` Ensembl does not know answered HTTP 400 on ``lookup/id``, and the fallback
written for an unknown *symbol* swallowed it: the gene background came out empty and the
run went on to a verdict, exit 0.  And a config whose transcripts a pinned release does not
have was refused for the gene background before the missing cDNA was named, which sent the
user to the config's gene when the release was the problem.
"""
import json
import time
import urllib.request

from isoform_dominance import cli
from test_ensembl_release import CURRENT, OLD, ONLY_NOW, SEQS, Ensembl, _http_error, _seq
from test_gene_choice import LEPR, Rest, _ids

UNKNOWN = "ENSG00000999999"


class RestWithoutGene(Rest):
    """lookup/id of an id Ensembl does not know answers 400, as the live server does."""

    def __call__(self, req, timeout=None, **_):
        if req.full_url.split("?")[0].endswith("/lookup/id/" + UNKNOWN):
            self.calls.append((req.get_method(), "/lookup/id/" + UNKNOWN, None))
            raise _http_error(req.full_url, 400, b'{"error":"ID not found"}')
        return super().__call__(req, timeout=timeout)


def _write(tmp_path, cfg):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(cfg))
    return str(p)


def test_a_gene_id_ensembl_does_not_know_stops(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(urllib.request, "urlopen", RestWithoutGene())
    tids = sorted(_ids(LEPR))
    cfg = _write(tmp_path, {"gene": "LEPR", "gene_id": UNKNOWN + ".3",
                            "groups": {"A": tids[:1], "B": tids[1:2]},
                            "primary_comparison": ["A", "B"]})
    assert cli.main(["identifiability", "--config", cfg]) == 1
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1, err
    assert "Ensembl release 116 has no gene %s" % UNKNOWN in err[0]


def test_an_unknown_symbol_still_leaves_the_background_empty(tmp_path, capsys, monkeypatch):
    # the fallback the gene_id path borrowed is right for a symbol, and stays
    class NoSymbol(Rest):
        def __call__(self, req, timeout=None, **_):
            if "/lookup/symbol/" in req.full_url:
                raise _http_error(req.full_url, 400)
            return super().__call__(req, timeout=timeout)
    monkeypatch.setattr(urllib.request, "urlopen", NoSymbol())
    tids = sorted(_ids(LEPR))
    cfg = _write(tmp_path, {"gene": "NOSUCHSYMBOL", "groups": {"A": tids[:1], "B": tids[1:2]},
                            "primary_comparison": ["A", "B"]})
    assert cli.main(["identifiability", "--config", cfg, "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["background"]["gene_transcripts"] == []


def test_transcripts_a_pinned_release_lacks_are_named_before_the_gene(tmp_path, capsys,
                                                                     monkeypatch):
    monkeypatch.setattr(time, "sleep", lambda s: None)
    second = "ENST%011d" % 151
    monkeypatch.setitem(SEQS[CURRENT], second, _seq(450, 98))
    monkeypatch.setattr(urllib.request, "urlopen", Ensembl())
    cfg = _write(tmp_path, {"gene": "FAKE", "groups": {"A": [ONLY_NOW], "B": [second]},
                            "primary_comparison": ["A", "B"], "ensembl_release": CURRENT})
    assert cli.main(["identifiability", "--config", cfg, "--ensembl-release", str(OLD)]) == 1
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1, err
    assert "Ensembl release %d has no cDNA for %s, %s" % (OLD, ONLY_NOW, second) in err[0]
