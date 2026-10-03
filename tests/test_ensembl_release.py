"""Pinning a run to an Ensembl release, with no network.

``rest.ensembl.org`` serves only the current release.  Earlier releases are served by
Ensembl's REST archive, behind an alias per release.  What the fake below reproduces was
observed against the live servers on 2026-09-24:

    e110.rest.ensembl.org   --301-->  jul2023.rest.ensembl.org    {"releases": [110]}
    e115.rest.ensembl.org   --301-->  sep2025.rest.ensembl.org    {"releases": [115]}
    e116.rest.ensembl.org   503; later the same day --301--> jun2026.rest.ensembl.org,
                            which answers 503 "Domain not existing": the current
                            release has no working archive alias
    e104.rest.ensembl.org   --301-->  may2021.rest.ensembl.org  --301-->
                            www.ensembl.org/help/articles/archives  (retired)
    rest.ensembl.org                                               {"releases": [116]}

    POST e110.rest.ensembl.org/sequence/id   ->  400 {"error":"ID '' not found"}
    POST jul2023.rest.ensembl.org/sequence/id ->  200, both sequences

The POST fails because ``urllib`` replays a POST that receives a 301, 302 or 303 as a
GET with no body.  So the alias can be used to *find* the archive host and never to fetch cDNA from
it, and the tests below assert on the host every request actually went to.
"""
import email.message
import io
import json
import random
import time
import urllib.request
from urllib.error import HTTPError

import pytest

from isoform_dominance import annotate, cli, ensembl, identifiability
from test_annotate_fixture import LEPR_LIKE

CURRENT, OLD = 116, 110
RETIRED_PAGE = "https://www.ensembl.org/help/articles/archives"
HOST = {CURRENT: "https://rest.ensembl.org", OLD: "https://jul2023.rest.ensembl.org"}
ALIAS = ensembl.ARCHIVE_ALIAS


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


GROUP_A = ["ENST%011d" % i for i in range(1, 9)]
GROUP_B = ["ENST%011d" % i for i in range(101, 104)]
BACKGROUND = ["ENST%011d" % i for i in range(201, 213)]
# The annotation grew between the two releases, and the same ids carry different
# sequence in each -- of a different length, because random sequences share no k-mer,
# so a class's unique k-mer count is only its window count and content alone would not
# show.  A run that read the wrong release therefore gives a different background and
# different counts, which is what lets the tests below fail.
BACKGROUND_AT = {CURRENT: BACKGROUND, OLD: BACKGROUND[:5]}
SEQS = {rel: {t: _seq(300 + 200 * (t in GROUP_B) + 20 * (rel == OLD), rel * 10_000 + i)
              for i, t in enumerate(GROUP_A + GROUP_B + BACKGROUND)}
        for rel in (CURRENT, OLD)}
# a transcript the current release has and the old one does not; a batch POST that names
# it gets the others back and nothing for it, as the live archive does
ONLY_NOW = "ENST%011d" % 150
SEQS[CURRENT][ONLY_NOW] = _seq(400, 99)
CONFIG = {"gene": "FAKE", "groups": {"A": GROUP_A, "B": GROUP_B},
          "primary_comparison": ["B", "A"]}

# LEPR as the old release annotated it: one short-class transcript not yet added
LEPR_OLD = dict(LEPR_LIKE, Transcript=[t for t in LEPR_LIKE["Transcript"]
                                       if t["id"] != "ENST00000616738"])


def _fake_lookup(rel):
    return {"id": "ENSG00000000001", "strand": 1,
            "Transcript": [{"id": "%s.%d" % (t, 1 + (rel == CURRENT))}
                           for t in GROUP_A + GROUP_B + BACKGROUND_AT[rel]]}


LOOKUPS = {CURRENT: {"FAKE": _fake_lookup(CURRENT), "LEPR": LEPR_LIKE},
           OLD: {"FAKE": _fake_lookup(OLD), "LEPR": LEPR_OLD}}


class _Resp(io.BytesIO):
    status = 200

    def __init__(self, payload, url):
        super().__init__(payload.encode() if isinstance(payload, str) else payload)
        self.headers = email.message.Message()
        self.url = url

    def geturl(self):
        return self.url


def _http_error(url, code, body=b""):
    return HTTPError(url, code, "status %d" % code, email.message.Message(), io.BytesIO(body))


class Ensembl:
    """rest.ensembl.org at ``CURRENT`` and the release-``OLD`` archive behind its alias.

    Every other alias answers 503, as ``e116`` and ``e117`` did live (``e116`` by way of
    a redirect to a host that answers 503).  ``calls`` logs
    ``(method, host, path)`` for every request, with the host it was *sent* to.
    """

    def __init__(self, alias_reports=OLD, alias_down=False, server_down=False,
                 alias_page=None, current_alias_page=None, alias_retired=False,
                 current_alias_to_rest=False):
        self.calls = []
        self.alias_reports = alias_reports
        self.alias_down = alias_down
        self.server_down = server_down
        self.alias_page = alias_page                  # a 200 that is not a listing
        self.current_alias_page = current_alias_page
        self.alias_retired = alias_retired            # redirects to the archives web page
        self.current_alias_to_rest = current_alias_to_rest   # e<current> -> rest.ensembl.org

    def __call__(self, req, timeout=None, **_):
        url, method = req.full_url, req.get_method()
        host = "/".join(url.split("/")[:3])
        path = url[len(host):]
        self.calls.append((method, host, path))
        if host == ALIAS % CURRENT and self.current_alias_to_rest:
            return self._serve(CURRENT, method, path, req, HOST[CURRENT] + path)
        if host == ALIAS % CURRENT and self.current_alias_page is not None:
            return _Resp(self.current_alias_page, url)
        if host == ALIAS % OLD:
            if self.alias_down:
                raise _http_error(url, 503)
            if self.alias_retired:           # as e104 and earlier answered on 2026-09-24
                return _Resp("<html>Ensembl archives</html>", RETIRED_PAGE)
            if self.alias_page is not None:
                return _Resp(self.alias_page, HOST[OLD] + path)
            if method == "POST":             # urllib's replay of a redirected POST
                raise _http_error(HOST[OLD] + path, 400, b'{"error":"ID \'\' not found"}')
            return self._serve(OLD, method, path, req, HOST[OLD] + path,
                               reports=self.alias_reports)
        if host == HOST[CURRENT] and self.server_down:
            raise _http_error(url, 503)
        for rel, h in HOST.items():
            if host == h:
                return self._serve(rel, method, path, req, url)
        raise _http_error(url, 503)

    @staticmethod
    def _serve(rel, method, path, req, final_url, reports=None):
        if path.startswith("/info/data"):
            return _Resp(json.dumps({"releases": [reports or rel]}), final_url)
        if path.startswith("/lookup/symbol/"):
            gene = path.split("/")[4].split("?")[0]
            if gene not in LOOKUPS[rel]:
                raise _http_error(final_url, 400)
            return _Resp(json.dumps(LOOKUPS[rel][gene]), final_url)
        if path.startswith("/xrefs/symbol/"):    # the one gene of the symbol, or none
            gene = path.split("/")[4].split("?")[0]
            return _Resp(json.dumps([{"type": "gene", "id": LOOKUPS[rel][gene]["id"]}]
                                    if gene in LOOKUPS[rel] else []), final_url)
        if method == "POST" and path.startswith("/sequence/id"):
            ids = [i.split(".")[0] for i in json.loads(req.data)["ids"]]
            return _Resp(json.dumps([{"query": i, "id": i, "seq": SEQS[rel][i]}
                                     for i in ids if i in SEQS[rel]]), final_url)
        raise AssertionError("unexpected request %s %s" % (method, path))

    def hosts(self, what=""):
        return [h for m, h, p in self.calls if p.startswith(what) or m == what]


@pytest.fixture
def sleeps(monkeypatch):
    waited = []
    monkeypatch.setattr(time, "sleep", waited.append)
    return waited


@pytest.fixture
def serve(monkeypatch, sleeps):
    def _serve(**kw):
        fake = Ensembl(**kw)
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        return fake
    return _serve


@pytest.fixture
def ens(serve):
    return serve()


@pytest.fixture
def no_network(monkeypatch):
    def refuse(req, *a, **k):
        raise AssertionError("network request %s" % req.full_url)
    monkeypatch.setattr(urllib.request, "urlopen", refuse)


def _offline(rel):
    return ({t: SEQS[rel][t] for t in GROUP_A + GROUP_B},
            {t: SEQS[rel][t] for t in BACKGROUND_AT[rel]})


# --------------------------------------------------------------------------- #
# resolving a release to a host
# --------------------------------------------------------------------------- #
def test_no_release_is_the_current_server_and_asks_nothing(no_network):
    assert ensembl.resolve_server(None) == ensembl.SERVER


def test_the_current_release_has_no_alias_and_is_rest_ensembl_org(ens, sleeps):
    assert ensembl.resolve_server(CURRENT) == HOST[CURRENT]
    assert ens.hosts() == [ALIAS % CURRENT, HOST[CURRENT]] and sleeps == []


def test_an_earlier_release_is_the_host_its_alias_redirects_to(ens):
    assert ensembl.resolve_server(OLD) == HOST[OLD]
    assert ens.calls == [("GET", ALIAS % OLD, "/info/data")]


def test_an_earlier_release_does_not_need_rest_ensembl_org(serve):
    serve(server_down=True)
    assert ensembl.resolve_server(OLD) == HOST[OLD]


def test_a_release_after_the_current_one_is_refused_without_retrying_the_archive(ens, sleeps):
    with pytest.raises(ensembl.ReleaseNotServed) as e:
        ensembl.resolve_server(CURRENT + 1)
    assert ens.hosts() == [ALIAS % (CURRENT + 1), HOST[CURRENT]] and sleeps == []
    # 116 is the last release on the REST API: "does not exist yet" would be false there
    assert "release 117 is not on the REST API" in str(e.value)
    assert "serves release 116, the last release" in str(e.value)
    assert "does not exist yet" not in str(e.value)


def test_a_release_after_a_current_one_below_the_last_rest_release_does_not_exist_yet(
        ens, monkeypatch):
    monkeypatch.setattr(ensembl, "REST_LAST_RELEASE", CURRENT + 5)
    with pytest.raises(ensembl.ReleaseNotServed, match="release 117 does not exist yet: "
                       "https://rest.ensembl.org serves 116, the current release"):
        ensembl.resolve_server(CURRENT + 1)


def test_a_retired_archive_is_named_as_retired_at_once(serve, sleeps):
    # a retired archive's alias redirects to the web page on archives: that is certain,
    # so no retry, no question to rest.ensembl.org, and no "down or retired" hedge
    fake = serve(alias_retired=True)
    with pytest.raises(ensembl.ReleaseNotServed) as e:
        ensembl.resolve_server(OLD)
    assert str(e.value) == ("the REST archive for Ensembl release 110 is retired: "
                            "https://e110.rest.ensembl.org redirects to %s, outside the "
                            "REST service" % RETIRED_PAGE)
    assert fake.hosts() == [ALIAS % OLD] and sleeps == []


def test_a_retired_archive_is_one_line_and_exit_1(serve, tmp_path, capsys):
    serve(alias_retired=True)
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(CONFIG))
    assert cli.main(["identify", "--config", str(cfg), "--ensembl-release", str(OLD)]) == 1
    err = capsys.readouterr().err
    assert err.startswith("Ensembl release not available: the REST archive for Ensembl "
                          "release 110 is retired") and len(err.strip().splitlines()) == 1


def test_an_archive_that_answers_for_another_release_is_refused(serve):
    serve(alias_reports=OLD + 1)
    with pytest.raises(ensembl.ReleaseNotServed, match="answered for release 111, not 110"):
        ensembl.resolve_server(OLD)


MAINTENANCE = b"<html><body>Ensembl is down for scheduled maintenance</body></html>"


def test_an_archive_answering_with_a_page_is_unreachable_not_a_config_error(serve):
    serve(alias_page=MAINTENANCE)
    with pytest.raises(ensembl.ReleaseNotServed, match="did not answer with a release listing"):
        ensembl.resolve_server(OLD)


def test_a_page_at_the_current_release_alias_still_resolves_to_rest_ensembl_org(serve):
    serve(current_alias_page=MAINTENANCE)
    assert ensembl.resolve_server(CURRENT) == HOST[CURRENT]


def test_an_alias_that_redirects_to_rest_ensembl_org_itself_is_not_retired(serve, sleeps):
    # rest.ensembl.org is inside the REST service, so the retirement test must accept it.
    # 116 is the release everyone will pin, its alias has no archive of its own, and
    # Ensembl could point it at the current server at any time; a host test that reads
    # only the ".rest.ensembl.org" suffix calls that a retired archive.
    fake = serve(current_alias_to_rest=True)
    assert ensembl.resolve_server(CURRENT) == HOST[CURRENT]
    assert fake.hosts() == [ALIAS % CURRENT] and sleeps == []


def test_an_unreachable_archive_is_reported_after_the_retries(serve, sleeps):
    serve(alias_down=True)
    with pytest.raises(ensembl.ReleaseNotServed, match="temporarily down or retired"):
        ensembl.resolve_server(OLD, retries=2, retry_wait=1.0)
    assert sleeps == [1.0, 2.0]


# --------------------------------------------------------------------------- #
# identifiability
# --------------------------------------------------------------------------- #
def test_cdna_is_never_posted_to_the_alias(ens):
    identifiability.analyze(CONFIG, ensembl_release=OLD)
    posts = ens.hosts("POST")
    assert posts and set(posts) == {HOST[OLD]}


def test_a_pinned_run_reads_sequence_and_background_from_that_release(ens):
    pinned = identifiability.analyze(CONFIG, ensembl_release=OLD)
    seqs, bg = _offline(OLD)
    offline = identifiability.analyze(CONFIG, sequences=seqs, background_sequences=bg)
    assert pinned["annotation"]["fetched_release"] == OLD
    assert pinned["background"]["gene_transcripts"] == sorted(BACKGROUND_AT[OLD])
    for g in ("A", "B"):
        assert pinned["groups"][g]["n_unique_kmers"] == offline["groups"][g]["n_unique_kmers"]
    # and the fixture can tell the releases apart, so the equality above means something
    now = identifiability.analyze(CONFIG)
    assert now["annotation"]["fetched_release"] == CURRENT
    assert now["background"]["n_background_transcripts"] == len(BACKGROUND)
    assert now["groups"]["A"]["n_unique_kmers"] != pinned["groups"]["A"]["n_unique_kmers"]


def test_every_request_of_a_pinned_run_after_resolving_goes_to_the_archive_host(ens):
    identifiability.analyze(CONFIG, ensembl_release=OLD)
    after = ens.calls[[c[1] for c in ens.calls].index(ALIAS % OLD) + 1:]
    assert after and {h for _, h, _ in after} == {HOST[OLD]}


def test_a_run_on_supplied_sequence_resolves_nothing(no_network):
    seqs, bg = _offline(OLD)
    res = identifiability.analyze(CONFIG, sequences=seqs, background_sequences=bg,
                                  ensembl_release=OLD)
    assert res["annotation"]["fetched_release"] is None


# --------------------------------------------------------------------------- #
# annotate
# --------------------------------------------------------------------------- #
def test_annotate_proposes_groups_from_the_release_asked_for(ens):
    old = annotate.build_config("LEPR", release=OLD)
    now = annotate.build_config("LEPR")
    assert (old["ensembl_release"], now["ensembl_release"]) == (OLD, CURRENT)
    assert "ENST00000616738" not in sum(old["groups"].values(), [])
    assert "ENST00000616738" in sum(now["groups"].values(), [])
    assert ens.hosts("/lookup/") == [HOST[OLD], HOST[CURRENT]]


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def test_cli_annotate_takes_the_release(ens, tmp_path):
    out = tmp_path / "cfg.json"
    assert cli.main(["annotate", "--gene", "LEPR", "--out", str(out),
                     "--ensembl-release", str(OLD)]) == cli.EXIT_OK
    assert json.loads(out.read_text())["ensembl_release"] == OLD


def test_cli_identify_takes_the_release(ens, tmp_path, capsys):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(dict(CONFIG, ensembl_release=OLD)))
    assert cli.main(["identify", "--config", str(cfg), "--ensembl-release", str(OLD),
                     "--json"]) == cli.EXIT_OK
    rep = json.loads(capsys.readouterr().out)
    assert rep["annotation"] == {"ensembl_release": OLD, "fetched_release": OLD,
                                 "file_release": None, "source": None}


def test_the_mismatch_note_names_the_flag_that_fetches_the_configs_release(ens, tmp_path,
                                                                           capsys):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(dict(CONFIG, ensembl_release=OLD)))
    cli.main(["identify", "--config", str(cfg)])
    assert "pass --ensembl-release %d" % OLD in capsys.readouterr().err


def test_a_deliberate_cross_release_run_says_what_it_did(ens, tmp_path, capsys):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(dict(CONFIG, ensembl_release=CURRENT)))
    cli.main(["identify", "--config", str(cfg), "--ensembl-release", str(OLD)])
    err = capsys.readouterr().err
    assert "fetched from Ensembl release %d, as requested" % OLD in err
    assert "proposed from release %d" % CURRENT in err


def test_a_transcript_the_pinned_release_lacks_is_named_with_the_release(ens, tmp_path,
                                                                        capsys):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(dict(CONFIG, ensembl_release=CURRENT,
                                   groups={"A": GROUP_A, "B": GROUP_B + [ONLY_NOW]})))
    assert cli.main(["identify", "--config", str(cfg), "--ensembl-release", str(OLD)]) == 1
    err = capsys.readouterr().err
    assert "Ensembl release %d has no cDNA for %s" % (OLD, ONLY_NOW) in err


@pytest.mark.parametrize("bad", ["0", "-5", "116.0", "latest"])
def test_the_release_must_be_a_positive_integer(no_network, tmp_path, capsys, bad):
    with pytest.raises(SystemExit) as e:
        cli.main(["annotate", "--gene", "LEPR", "--out", str(tmp_path / "o.json"),
                  "--ensembl-release", bad])
    assert e.value.code == 2 and "positive release number" in capsys.readouterr().err


@pytest.mark.parametrize("cmd", ["annotate", "identify"])
def test_an_unavailable_release_is_one_line_and_exit_1_not_a_traceback(ens, tmp_path,
                                                                        capsys, cmd):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(CONFIG))
    argv = (["annotate", "--gene", "LEPR", "--out", str(tmp_path / "o.json")]
            if cmd == "annotate" else ["identify", "--config", str(cfg)])
    assert cli.main(argv + ["--ensembl-release", str(CURRENT + 1)]) == 1
    err = capsys.readouterr().err
    assert err.startswith("Ensembl release not available:") and "not on the REST API" in err
