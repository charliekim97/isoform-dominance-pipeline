"""The Ensembl release is recorded, reported, and never claimed for sequence it did not
produce.

The verdict is a function of the annotation release -- between GENCODE v44 and Ensembl
116, 6 of 36 verdicts in a 49-gene panel moved -- so a config that does not say which
release its groups came from gives a verdict nobody can reproduce.  ``rest.ensembl.org``
serves only the current release, so the release a config was annotated against and the
release ``identifiability`` later fetches sequence from are different facts, and are
reported separately.
"""
import json
import urllib.request

import pytest

from isoform_dominance import annotate, cli, identifiability
from test_annotate_fixture import LEPR_LIKE
from test_ensembl_http import BACKGROUND, CONFIG, GROUP_A, GROUP_B, SEQS, FakeEnsembl

OFFLINE = {t: SEQS[t] for t in GROUP_A + GROUP_B}


@pytest.fixture
def serve(monkeypatch):
    def _serve(**kw):
        fake = FakeEnsembl(**kw)
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        return fake
    return _serve


@pytest.fixture
def no_network(monkeypatch):
    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)


# --------------------------------------------------------------------------- #
# annotate
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("release", [115, 116])
def test_annotate_records_the_release_it_fetched_from(serve, release):
    fake = serve(lookup=LEPR_LIKE, release=release)
    cfg = annotate.build_config("LEPR")
    assert cfg["ensembl_release"] == release
    assert "/info/data" in [c[1] for c in fake.calls]


def test_annotate_no_longer_writes_the_free_text_reference(serve):
    serve(lookup=LEPR_LIKE)
    assert "reference" not in annotate.build_config("LEPR")


# --------------------------------------------------------------------------- #
# identifiability
# --------------------------------------------------------------------------- #
def test_the_report_carries_the_config_release(no_network):
    res = identifiability.analyze(dict(CONFIG, ensembl_release=116), sequences=OFFLINE)
    assert res["annotation"] == {"ensembl_release": 116, "fetched_release": None}


def test_an_old_config_reports_no_release_and_still_runs(no_network):
    old = dict(CONFIG, reference="Ensembl REST (live annotation)")
    res = identifiability.analyze(old, sequences=OFFLINE)
    assert res["annotation"] == {"ensembl_release": None, "fetched_release": None}


def test_a_run_that_fetches_sequence_records_the_release_it_fetched_from(serve):
    """Annotated against 116, run against a server now at 117: both are reported, and
    the config's release is not passed off as the release of the sequence used."""
    serve(release=117)
    res = identifiability.analyze(dict(CONFIG, ensembl_release=116))
    assert res["annotation"] == {"ensembl_release": 116, "fetched_release": 117}


def test_a_pinned_run_makes_no_request_at_all(no_network, tmp_path):
    """--sequences and --background-fasta from one release is how a verdict is pinned.
    --background-fasta alone is not: the configured transcripts' cDNA and the gene
    background are still fetched live.  Asking the server for its release here would be
    a request, and a claim, about sequence this run did not use."""
    fa = tmp_path / "bg.fa"
    fa.write_text(">DECOY\n%s\n" % SEQS[BACKGROUND[0]])
    res = identifiability.analyze(dict(CONFIG, ensembl_release=116),
                                  sequences=OFFLINE, background_fasta=str(fa))
    assert res["annotation"]["fetched_release"] is None


# --------------------------------------------------------------------------- #
# the CLI
# --------------------------------------------------------------------------- #
def _cli_case(tmp_path, **extra):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(dict(CONFIG, **extra)))
    sq = tmp_path / "seqs.json"
    sq.write_text(json.dumps(OFFLINE))
    return ["identify", "--config", str(cfg), "--sequences", str(sq)]


def test_the_header_names_the_release(no_network, tmp_path, capsys):
    assert cli.main(_cli_case(tmp_path, ensembl_release=116)) == cli.EXIT_OK
    out = capsys.readouterr()
    assert "annotation: Ensembl release 116" in out.out
    assert "not reproducible" not in out.err


def test_a_config_without_a_release_says_so_and_does_not_fail(no_network, tmp_path, capsys):
    argv = _cli_case(tmp_path, reference="Ensembl REST (live annotation)")
    assert cli.main(argv) == cli.EXIT_OK
    out = capsys.readouterr()
    assert "annotation: release not recorded in the config" in out.out
    assert len([ln for ln in out.err.splitlines() if "not reproducible" in ln]) == 1


def test_a_release_mismatch_is_reported(serve, tmp_path, capsys):
    serve(release=117)
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(dict(CONFIG, ensembl_release=116)))
    assert cli.main(["identify", "--config", str(cfg)]) == cli.EXIT_OK
    out = capsys.readouterr()
    assert "annotation: Ensembl release 116; sequence fetched from release 117" in out.out
    assert "116" in out.err and "117" in out.err


def test_the_json_report_carries_the_release(no_network, tmp_path, capsys):
    cli.main(_cli_case(tmp_path, ensembl_release=116) + ["--json"])
    assert json.loads(capsys.readouterr().out)["annotation"]["ensembl_release"] == 116
