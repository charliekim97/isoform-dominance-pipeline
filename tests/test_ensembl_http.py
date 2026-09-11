"""Ensembl access with no network: batched cDNA, bounded retries, no silent partials.

The HTTP layer is replaced at ``urllib.request.urlopen`` and the clock at
``time.sleep``, so these tests exercise the real request construction, retry loop and
response parsing, and assert on the requests actually made and the waits actually taken.
"""
import email.message
import io
import json
import random
import time
import urllib.request
from urllib.error import HTTPError

import pytest

from isoform_dominance import annotate, cli, identifiability
from test_annotate_fixture import LEPR_LIKE

SERVER = "https://rest.ensembl.org"


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


# 60 + 5 configured transcripts and 40 more of the same gene, as background.  The three
# sets have different lengths, so handing one transcript another's sequence changes the
# k-mer counts -- which is what lets the equivalence test below fail.
GROUP_A = ["ENST%011d" % i for i in range(1, 61)]
GROUP_B = ["ENST%011d" % i for i in range(101, 106)]
BACKGROUND = ["ENST%011d" % i for i in range(201, 241)]
ALL_IDS = GROUP_A + GROUP_B + BACKGROUND
SEQS = {t: _seq({0: 300, 1: 500, 2: 400}[(t in GROUP_B) + 2 * (t in BACKGROUND)], i)
        for i, t in enumerate(ALL_IDS)}
CONFIG = {"gene": "FAKE", "groups": {"A": GROUP_A, "B": GROUP_B},
          "primary_comparison": ["B", "A"]}
GENE_LOOKUP = {"id": "ENSG00000000001", "strand": 1,
               "Transcript": [{"id": t + ".1"} for t in ALL_IDS]}


class _Resp(io.BytesIO):
    status = 200

    def __init__(self, payload):
        super().__init__(payload.encode() if isinstance(payload, str) else payload)
        self.headers = email.message.Message()


def _http_error(url, code, retry_after=None):
    hdrs = email.message.Message()
    if retry_after is not None:
        hdrs["Retry-After"] = str(retry_after)
    return HTTPError(url, code, "status %d" % code, hdrs, io.BytesIO(b""))


class FakeEnsembl:
    """Serves lookups and cDNA, one id per GET or up to 50 per POST, and logs each call."""

    def __init__(self, fail=None, lookup=GENE_LOOKUP):
        self.calls = []
        self.fail = fail or (lambda call: None)
        self.lookup = lookup

    def __call__(self, req, timeout=None, **_):
        path = req.full_url[len(SERVER):]
        body = json.loads(req.data) if req.data else None
        call = (req.get_method(), path, body)
        self.calls.append(call)
        err = self.fail(call)
        if err is not None:
            raise err
        if path.startswith("/lookup/symbol/"):
            return _Resp(json.dumps(self.lookup))
        if call[0] == "POST" and path.startswith("/sequence/id"):
            ids = [i.split(".")[0] for i in body["ids"]]
            # reversed, so a caller that pairs results with ids by position is caught
            return _Resp(json.dumps([{"query": i, "id": i, "seq": SEQS[i], "molecule": "dna"}
                                     for i in reversed(ids)]))
        if call[0] == "GET" and path.startswith("/sequence/id/"):
            return _Resp(SEQS[path.split("/")[3].split("?")[0]])
        raise AssertionError("unexpected request %r" % (call[:2],))

    def sequence_calls(self):
        return [c for c in self.calls if c[1].startswith("/sequence/id")]


def _ids_in(call):
    method, path, body = call
    return body["ids"] if method == "POST" else [path.split("/")[3].split("?")[0]]


def fail_first(n, make, on="/sequence/id"):
    left = [n]

    def fail(call):
        if call[1].startswith(on) and left[0] > 0:
            left[0] -= 1
            return make(SERVER + call[1])
        return None
    return fail


@pytest.fixture
def sleeps(monkeypatch):
    waited = []
    monkeypatch.setattr(time, "sleep", waited.append)
    return waited


@pytest.fixture
def serve(monkeypatch):
    def _serve(fail=None, lookup=GENE_LOOKUP):
        fake = FakeEnsembl(fail, lookup)
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        return fake
    return _serve


# --------------------------------------------------------------------------- #
# batching
# --------------------------------------------------------------------------- #
def test_cdna_is_fetched_in_posts_of_at_most_50_ids(serve, sleeps):
    fake = serve()
    identifiability.analyze(CONFIG)
    seq = fake.sequence_calls()
    assert [m for m, _, _ in seq] == ["POST"] * len(seq), [c[:2] for c in seq[:3]]
    assert all("type=cdna" in p for _, p, _ in seq)
    assert all(len(b["ids"]) <= 50 for _, _, b in seq)
    assert len(seq) == 3                         # 65 configured + 40 background ids
    assert sorted(i for c in seq for i in _ids_in(c)) == sorted(ALL_IDS)
    assert sleeps == []


def test_batched_cdna_gives_the_answer_supplied_sequences_give(serve, sleeps):
    serve()
    online = identifiability.analyze(CONFIG)
    offline = identifiability.analyze(
        CONFIG, sequences={t: SEQS[t] for t in GROUP_A + GROUP_B},
        background_sequences={t: SEQS[t] for t in BACKGROUND})
    assert online["background"]["n_background_transcripts"] == len(BACKGROUND)
    for g in ("A", "B"):
        assert online["groups"][g]["n_unique_kmers"] == offline["groups"][g]["n_unique_kmers"]


# --------------------------------------------------------------------------- #
# retries
# --------------------------------------------------------------------------- #
def test_a_503_is_retried_with_exponential_backoff(serve, sleeps):
    serve(fail_first(2, lambda u: _http_error(u, 503)))
    res = identifiability.analyze(CONFIG)
    assert res["background"]["n_background_transcripts"] == len(BACKGROUND)
    assert sleeps == [1.0, 2.0]                  # the documented default: 1 s, doubling


def test_a_read_timeout_is_retried(serve, sleeps):
    serve(fail_first(1, lambda u: TimeoutError("The read operation timed out")))
    identifiability.analyze(CONFIG)
    assert sleeps == [1.0]


def test_429_waits_as_long_as_retry_after_says(serve, sleeps):
    serve(fail_first(1, lambda u: _http_error(u, 429, retry_after=7)))
    identifiability.analyze(CONFIG)
    assert sleeps == [7.0]


def test_429_without_retry_after_falls_back_to_the_backoff(serve, sleeps):
    serve(fail_first(1, lambda u: _http_error(u, 429)))
    identifiability.analyze(CONFIG)
    assert sleeps == [1.0]


def test_retry_count_and_wait_are_parameters_and_the_error_surfaces(serve, sleeps):
    fake = serve(fail_first(99, lambda u: _http_error(u, 503)))
    with pytest.raises(HTTPError) as e:
        identifiability.analyze(CONFIG, retries=2, retry_wait=0.25)
    assert e.value.code == 503
    assert len(fake.sequence_calls()) == 3       # the first try and two retries
    assert sleeps == [0.25, 0.5]


def test_a_client_error_is_not_retried(serve, sleeps):
    fake = serve(fail_first(99, lambda u: _http_error(u, 400)))
    with pytest.raises(HTTPError):
        identifiability.analyze(CONFIG)
    assert len(fake.sequence_calls()) == 1
    assert sleeps == []


# --------------------------------------------------------------------------- #
# the background is complete or the run fails
# --------------------------------------------------------------------------- #
def test_a_background_that_cannot_be_fetched_is_an_error_not_a_smaller_background(
        serve, sleeps):
    """The per-transcript path swallowed any exception part-way through the background
    and carried on with whatever had arrived, so a flaky network silently produced a
    different answer.  Here the first background transcript arrives and the rest fail."""
    rest = set(BACKGROUND[1:])
    serve(lambda call: (_http_error(SERVER + call[1], 503)
                        if call[1].startswith("/sequence/id")
                        and rest & {i.split(".")[0] for i in _ids_in(call)} else None))
    with pytest.raises(HTTPError):
        identifiability.analyze(CONFIG)


def test_an_unknown_gene_symbol_still_means_no_gene_background(serve, sleeps):
    serve(fail_first(1, lambda u: _http_error(u, 400), on="/lookup/"))
    res = identifiability.analyze(CONFIG)
    assert res["background"]["n_background_transcripts"] == 0
    assert sleeps == []


# --------------------------------------------------------------------------- #
# annotate, and the CLI
# --------------------------------------------------------------------------- #
def test_annotate_lookup_is_retried(serve, sleeps):
    serve(fail_first(1, lambda u: _http_error(u, 503), on="/lookup/"), lookup=LEPR_LIKE)
    cfg = annotate.build_config("LEPR")
    assert "ENST00000349533" in {t for ids in cfg["groups"].values() for t in ids}
    assert sleeps == [1.0]


def test_cli_exposes_the_retry_parameters(serve, sleeps, tmp_path):
    serve(fail_first(99, lambda u: _http_error(u, 503)))
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps(CONFIG))
    rc = cli.main(["identifiability", "--config", str(cfg),
                   "--retries", "1", "--retry-wait", "0.5"])
    assert rc == 1                               # reported as a network failure
    assert sleeps == [0.5]

    serve(fail_first(99, lambda u: _http_error(u, 503), on="/lookup/"), lookup=LEPR_LIKE)
    rc = cli.main(["annotate", "--gene", "LEPR", "--out", str(tmp_path / "a.json"),
                   "--retries", "0"])
    assert rc == 1
    assert sleeps == [0.5]                       # --retries 0: nothing further waited
