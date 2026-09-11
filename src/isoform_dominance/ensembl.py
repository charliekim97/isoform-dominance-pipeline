"""Ensembl REST access: batched sequence retrieval and bounded retries.

Every network request the package makes goes through :func:`request`.  A run of
``identifiability`` needs the cDNA of every configured transcript and of the gene's
other transcripts -- commonly 20 to 60 of them -- and fetching them one GET at a time
meant that one transient failure among dozens of requests killed the run.  Sequences
are fetched with ``POST /sequence/id`` instead, :data:`MAX_POST_IDS` ids per request,
and every request is retried on the failures worth retrying.

Retried
    HTTP 429 (waiting for the server's ``Retry-After`` when it sends one), HTTP 500,
    502, 503 and 504, connection errors, and read timeouts.
Not retried
    Any other HTTP status.  A 400 for an unknown id or gene symbol will not improve by
    asking again.
Wait
    ``retry_wait * 2**i`` seconds before retry ``i`` (counting from 0), unless a 429
    named its own wait.

Defaults: :data:`DEFAULT_RETRIES` = 5 retries after the first attempt,
:data:`DEFAULT_RETRY_WAIT` = 1.0 s, so up to 1 + 2 + 4 + 8 + 16 = 31 s of backoff per
request before the last error is raised; :data:`DEFAULT_TIMEOUT` = 30 s per attempt.
"""
import json
import time
import urllib.request
from urllib.error import HTTPError, URLError

SERVER = "https://rest.ensembl.org"

#: Ensembl's limit on ids per ``POST /sequence/id`` request.
MAX_POST_IDS = 50

#: Retries after the first attempt.
DEFAULT_RETRIES = 5

#: Seconds before the first retry; each later retry doubles it.
DEFAULT_RETRY_WAIT = 1.0

#: Seconds an attempt may wait for a response before it counts as failed.
DEFAULT_TIMEOUT = 30

#: HTTP statuses that are retried.
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})


def _retry_after(err):
    """Seconds a response's ``Retry-After`` header asks for, or None."""
    value = err.headers.get("Retry-After") if err.headers is not None else None
    try:
        return max(0.0, float(value))
    except (TypeError, ValueError):
        return None


def request(path, body=None, *, accept="application/json", timeout=DEFAULT_TIMEOUT,
            retries=DEFAULT_RETRIES, retry_wait=DEFAULT_RETRY_WAIT):
    """GET ``path``, or POST ``body`` to it as JSON, and return the response bytes.

    The last error is re-raised once ``retries`` retries have been spent.
    """
    if retries < 0:
        raise ValueError("retries must be >= 0, got %r" % (retries,))
    data = None if body is None else json.dumps(body).encode()
    headers = {"Content-Type": "application/json", "Accept": accept}
    for attempt in range(retries + 1):
        req = urllib.request.Request(SERVER + path, data=data, headers=headers,
                                     method="GET" if data is None else "POST")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except HTTPError as e:                  # before URLError: it is a subclass
            if e.code not in RETRY_STATUS or attempt == retries:
                raise
            wait = _retry_after(e) if e.code == 429 else None
        except (URLError, TimeoutError):
            if attempt == retries:
                raise
            wait = None
        time.sleep(retry_wait * 2 ** attempt if wait is None else wait)
    raise AssertionError("unreachable")


def get_json(path, **retry):
    """GET ``path`` and decode the JSON response; ``retry`` goes to :func:`request`."""
    return json.loads(request(path, **retry))


def fetch_cdna_batch(ids, **retry):
    """``{id: cdna}`` for Ensembl transcript ids, :data:`MAX_POST_IDS` per request.

    Versions are stripped from the ids.  Results are matched to ids by the ``query``
    field Ensembl echoes back (``id`` if it is absent), never by position.  An id for
    which Ensembl returns nothing is absent from the result; the caller decides whether
    that is an error.  ``retry`` goes to :func:`request`.
    """
    want = list(dict.fromkeys(t.split(".")[0] for t in ids))
    out = {}
    for i in range(0, len(want), MAX_POST_IDS):
        chunk = want[i:i + MAX_POST_IDS]
        for item in json.loads(request("/sequence/id?type=cdna", {"ids": chunk}, **retry)):
            out[(item.get("query") or item["id"]).split(".")[0]] = item["seq"]
    return out
