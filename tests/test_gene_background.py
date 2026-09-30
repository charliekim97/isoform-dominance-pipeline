"""The gene background is the config's gene, and saved inputs say which gene that was.

`identifiability` fetched a config's gene background by symbol.  For CD99 the symbol
gives the chrY copy, whose sequence is the chrX copy's, so a config of the chrX gene was
judged against its own sequence and had no unique k-mer: `not_identifiable`, against
`identifiable` with the same sequence and no gene background.  A config that records
``gene_id`` -- `annotate` now writes it -- is looked up by id instead, and one without it
is refused when the symbol's gene holds none of the transcripts it names.
"""
import json
import urllib.request

import pytest

from isoform_dominance import cli, identifiability, io
from test_gene_choice import CD99_X, CD99_Y, LEPR, Rest, _ids


@pytest.fixture
def rest(monkeypatch):
    def _serve(**kw):
        fake = Rest(**kw)
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        return fake
    return _serve


# --------------------------------------------------------------------------- #
# identifiability: the gene background is the config's gene
# --------------------------------------------------------------------------- #
def _cfg(gene, **extra):
    tids = sorted(_ids(gene))
    return dict({"gene": gene["display_name"],
                 "groups": {"A": tids[:1], "B": tids[1:2]}, "primary_comparison": ["A", "B"]},
                **extra)


def test_the_gene_background_is_fetched_by_the_config_gene_id(rest):
    fake = rest()
    res = identifiability.analyze(_cfg(CD99_X, gene_id=CD99_X["id"]))
    assert ("GET", "/lookup/id/%s" % CD99_X["id"]) in fake.paths()
    assert not any(p.startswith("/lookup/symbol") for _, p in fake.paths())
    assert set(res["background"]["gene_transcripts"]) == _ids(CD99_X) - set(sorted(_ids(CD99_X))[:2])
    # the chrX transcripts judged against chrX's own others: distinguishable
    assert all(res["groups"][g]["n_unique_kmers"] > 0 for g in ("A", "B"))


def test_a_symbol_background_that_holds_none_of_the_configured_transcripts_is_refused(
        rest, tmp_path, capsys):
    rest()
    with pytest.raises(ValueError, match="gene_id"):
        identifiability.analyze(_cfg(CD99_X))              # no gene_id; symbol gives chrY
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(_cfg(CD99_X)))
    assert cli.main(["identifiability", "--config", str(cfg)]) == 1
    assert CD99_Y["id"] in capsys.readouterr().err


def test_a_config_without_gene_id_still_uses_the_symbol(rest):
    fake = rest()
    identifiability.analyze(_cfg(LEPR))
    assert ("GET", "/lookup/symbol/homo_sapiens/LEPR") in fake.paths()


# --------------------------------------------------------------------------- #
# saved inputs carry the gene id
# --------------------------------------------------------------------------- #
def test_saved_inputs_record_the_gene_id_and_a_rerun_for_another_is_refused(
        rest, tmp_path, capsys):
    rest()
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(_cfg(CD99_X, gene_id=CD99_X["id"])))
    saved = tmp_path / "in.json"
    assert cli.main(["identifiability", "--config", str(cfg),
                     "--save-inputs", str(saved)]) != 1
    assert json.load(open(saved))["gene_id"] == CD99_X["id"]

    other = tmp_path / "other.json"
    other.write_text(json.dumps(_cfg(CD99_X, gene_id=CD99_Y["id"])))
    capsys.readouterr()
    assert cli.main(["identifiability", "--config", str(other), "--inputs", str(saved)]) == 1
    err = capsys.readouterr().err
    assert CD99_X["id"] in err and CD99_Y["id"] in err


def test_a_saved_gene_id_that_is_not_a_string_is_refused(rest, tmp_path):
    rest()
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(_cfg(LEPR)))
    saved = tmp_path / "in.json"
    cli.main(["identifiability", "--config", str(cfg), "--save-inputs", str(saved)])
    doc = json.load(open(saved))
    doc["gene_id"] = 7
    saved.write_text(json.dumps(doc))
    with pytest.raises(io.InputError, match="gene_id"):
        io.load_inputs(str(saved))
