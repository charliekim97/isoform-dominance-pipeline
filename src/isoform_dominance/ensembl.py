"""Ensembl REST access: batched sequence retrieval and bounded retries.

Every network request the package makes goes through :func:`request`.  A run of
``identifiability`` needs the cDNA of every configured transcript and of the gene's
other transcripts -- commonly 20 to 60 of them -- and fetching them one GET at a time
meant that one transient failure among dozens of requests killed the run.  Sequences
are fetched with ``POST /sequence/id`` instead, :data:`MAX_POST_IDS` ids per request,
and every request is retried on the failures worth retrying.

Retried
    HTTP 429 (waiting for the server's ``Retry-After`` when it sends one), HTTP 500,
    502, 503 and 504, connection errors, read timeouts, and the failures below HTTP that
    ``urllib`` does not wrap in ``URLError`` (:data:`TRANSIENT`): a connection closed
    before the status line, a garbled status line, a body cut short, a reset mid-body.
    A JSON request whose 200 body is not JSON is retried too.
Not retried
    Any other HTTP status.  A 400 for an unknown id or gene symbol will not improve by
    asking again.  Nor will a URL ``http.client`` refuses to send (``InvalidURL``).
Wait
    ``retry_wait * 2**i`` seconds before retry ``i`` (counting from 0), unless a 429
    named its own wait.

Defaults: :data:`DEFAULT_RETRIES` = 5 retries after the first attempt,
:data:`DEFAULT_RETRY_WAIT` = 1.0 s, so up to 1 + 2 + 4 + 8 + 16 = 31 s of backoff per
request before the last error is raised; :data:`DEFAULT_TIMEOUT` = 30 s per attempt.

Releases
    ``rest.ensembl.org`` serves only the current release, and release 116 is the last
    one it will serve (:data:`REST_LAST_RELEASE`).  Earlier releases are served
    by Ensembl's REST archive, reached through an alias -- ``e110.rest.ensembl.org``
    for release 110 -- that answers with an HTTP 301 to a date-named host
    (``jul2023.rest.ensembl.org``).  :func:`resolve_server` follows that redirect once,
    with a GET, and every later request goes to the host it lands on.  The alias is
    never used directly for the batched cDNA fetch: ``urllib`` replays a POST that
    receives a 301, 302 or 303 as a GET with no body, and Ensembl answers that with
    ``400 {"error":"ID '' not found"}``.  The release a host reports is checked against
    the one asked for before anything is fetched from it.  An archive Ensembl has retired
    redirects off the REST service instead, to its web page on archives
    (``www.ensembl.org/help/articles/archives``).  On 2026-09-24 releases 90 to 104 did
    so, except 94; 94 and 75 to 89 answered 503; 105 onward answered, 106 and 111
    intermittently.  A redirect that ends anywhere else off the REST service is retried
    like an outage: only that page is taken for a retirement.
"""
import http.client
import json
import ssl
import time
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError

SERVER = "https://rest.ensembl.org"

#: Alias for the REST archive of one release; it redirects to a date-named host.
ARCHIVE_ALIAS = "https://e%d.rest.ensembl.org"

#: The last Ensembl release published on the REST API.  Ensembl 116 (June 2026) is the
#: final release of the legacy platform; its REST API "remains available for e116 for
#: long term use" and has "no plans to port over" to the new platform, which publishes
#: later releases through GraphQL and refget only (ensembl.info, 2026-07-21).  Used only
#: to word the error for a later release; which release is current is always asked.
REST_LAST_RELEASE = 116

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

#: Failures worth another attempt, and what a caller should catch as "the network failed".
#: ``urllib`` wraps errors raised while *sending* a request in ``URLError``, but not those
#: raised while reading the answer, so these arrive raw -- measured against a local socket
#: server: a connection closed before any response is ``RemoteDisconnected``, a garbled
#: status line ``BadStatusLine``, a body shorter than its ``Content-Length``
#: ``IncompleteRead``, a reset mid-body ``ConnectionResetError``.  Before 2.4 they were
#: neither retried nor reported: one attempt, then a traceback.  ``HTTPError`` is a
#: ``URLError`` and so is included, but it is retried only for :data:`RETRY_STATUS`.
TRANSIENT = (URLError, TimeoutError, ConnectionError, http.client.HTTPException,
             ssl.SSLError)


class BadResponse(URLError):
    """A 200 whose body is not what the request asked for, such as HTML where JSON was due.

    A ``URLError``, so it is retried and reported like any other network failure rather
    than surfacing as a JSON decoding error that reads like a config problem.
    """


def _retry_after(err):
    """Seconds a response's ``Retry-After`` header asks for, or None."""
    value = err.headers.get("Retry-After") if err.headers is not None else None
    try:
        return max(0.0, float(value))
    except (TypeError, ValueError):
        return None


def _with_retries(url, data, headers, read, *, timeout, retries, retry_wait):
    """Open ``url`` under the retry policy above and return ``read(response)``."""
    if retries < 0:
        raise ValueError("retries must be >= 0, got %r" % (retries,))
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, data=data, headers=headers,
                                     method="GET" if data is None else "POST")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return read(resp)
        except HTTPError as e:                  # before URLError: it is a subclass
            if e.code not in RETRY_STATUS or attempt == retries:
                raise
            wait = _retry_after(e) if e.code == 429 else None
        except http.client.InvalidURL:          # an HTTPException, but asking again cannot fix it
            raise
        except TRANSIENT:
            if attempt == retries:
                raise
            wait = None
        time.sleep(retry_wait * 2 ** attempt if wait is None else wait)
    raise AssertionError("unreachable")


def _decode_json(raw, url):
    try:
        return json.loads(raw)
    except ValueError as e:
        head = raw[:80].decode("utf-8", "replace") if isinstance(raw, bytes) else str(raw)[:80]
        raise BadResponse("%s answered with something that is not JSON: %r"
                          % (url, head)) from e


def request(path, body=None, *, server=None, accept="application/json",
            timeout=DEFAULT_TIMEOUT, retries=DEFAULT_RETRIES, retry_wait=DEFAULT_RETRY_WAIT):
    """GET ``path``, or POST ``body`` to it as JSON, and return the response bytes.

    ``server`` is the base URL, :data:`SERVER` when None; for an earlier release pass
    what :func:`resolve_server` returns, never the ``e<N>`` alias itself.  The last
    error is re-raised once ``retries`` retries have been spent.
    """
    return _request(path, body, lambda r: r.read(), server=server, accept=accept,
                    timeout=timeout, retries=retries, retry_wait=retry_wait)


def _request(path, body, read, *, server=None, accept="application/json",
             timeout=DEFAULT_TIMEOUT, retries=DEFAULT_RETRIES, retry_wait=DEFAULT_RETRY_WAIT):
    data = None if body is None else json.dumps(body).encode()
    headers = {"Content-Type": "application/json", "Accept": accept}
    return _with_retries((server or SERVER) + path, data, headers, read,
                         timeout=timeout, retries=retries, retry_wait=retry_wait)


def request_json(path, body=None, **retry):
    """As :func:`request`, decoding the JSON answer inside the retry loop.

    A 200 that is not JSON raises :class:`BadResponse` and is retried like any other
    transient failure.
    """
    url = (retry.get("server") or SERVER) + path
    return _request(path, body, lambda r: _decode_json(r.read(), url), **retry)


def get_json(path, **retry):
    """GET ``path`` and decode the JSON response; ``retry`` goes to :func:`request`."""
    return request_json(path, **retry)


def release_number(info):
    """The release in an ``/info/data`` response, which lists the releases a server holds.

    ``rest.ensembl.org`` holds one.  Anything else is raised rather than guessed at,
    because the number is recorded as the release a verdict was computed on.
    """
    releases = info.get("releases") or []
    if len(releases) != 1:
        raise ValueError("Ensembl /info/data listed releases %r; expected exactly one"
                         % (releases,))
    return int(releases[0])


def fetch_release(**retry):
    """The Ensembl release the server is serving; ``retry`` goes to :func:`request`."""
    return release_number(get_json("/info/data", **retry))


class ReleaseNotServed(LookupError):
    """No REST server answers for the Ensembl release that was asked for."""


class _NotAListing(Exception):
    """An archive answered, but not with a one-release ``/info/data`` listing."""


#: The host of :data:`SERVER`, which is inside the REST service but is not a subdomain of
#: itself: an archive alias that lands here has been pointed at the current server, not
#: retired.  Release 116 is the release everyone will pin and its alias has no archive of
#: its own, so this is the redirect to expect.
_REST_HOST = urllib.parse.urlsplit(SERVER).hostname


def _in_rest_service(host):
    """Is ``host`` part of the REST service, rather than the web site an archive retires to?"""
    return host == _REST_HOST or (host or "").endswith("." + _REST_HOST)


#: Where Ensembl sends a retired REST archive: its web page on archives.
RETIRED_HOSTS = frozenset({"www.ensembl.org", "ensembl.org"})
RETIRED_PATH = "/help/articles/archives"


def _retired(url):
    """Is ``url`` the page a retired REST archive redirects to?"""
    parts = urllib.parse.urlsplit(url)
    return parts.hostname in RETIRED_HOSTS and parts.path.startswith(RETIRED_PATH)


class _OffService(URLError):
    """A redirect that ended outside the REST service, but not on the page on archives: a
    maintenance page, say.  A ``URLError``, so it is retried like an outage."""


def _archive(alias, release, **net):
    """The host ``alias`` redirects to, provided it reports exactly ``release``.

    A 200 that is not a one-release listing -- a maintenance page, say -- raises
    :class:`_NotAListing`, which :func:`resolve_server` treats as an archive it could not
    reach rather than as a config error.  A redirect to Ensembl's page on archives
    (:data:`RETIRED_HOSTS`, a path under :data:`RETIRED_PATH`) is a retired archive, and
    raises :class:`ReleaseNotServed` at once.  One that ends anywhere else outside the REST
    service is retried under ``net``'s policy, and then raises :class:`_OffService`.
    """
    def _read(r):
        final = r.geturl()
        if _retired(final):
            # how Ensembl retires an archive (e100, e103, e104 on 2026-10-02, and 90-104
            # but 94 on 2026-09-24): certain, so not retried
            raise ReleaseNotServed("the REST archive for Ensembl release %d is retired: %s "
                                   "redirects to %s, outside the REST service"
                                   % (release, alias, final))
        if not _in_rest_service(urllib.parse.urlsplit(final).hostname):
            # a maintenance page for an hour is not a retirement
            raise _OffService("%s redirects to %s, outside the REST service"
                              % (alias, final))
        return r.read(), final

    raw, final = _with_retries(alias + "/info/data", None, {"Accept": "application/json"},
                               _read, **net)
    base = final.split("/info/data", 1)[0]
    try:
        got = release_number(json.loads(raw))
    except (ValueError, AttributeError, TypeError) as e:
        raise _NotAListing("%s did not answer with a release listing (%s)" % (base, e)) from e
    if got != release:
        raise ReleaseNotServed("%s answered for release %d, not %d" % (base, got, release))
    return base


def resolve_server(release=None, *, timeout=DEFAULT_TIMEOUT, retries=DEFAULT_RETRIES,
                   retry_wait=DEFAULT_RETRY_WAIT):
    """The REST base URL that serves Ensembl ``release``; :data:`SERVER` when None.

    An earlier release is reached through :data:`ARCHIVE_ALIAS`, resolved here with a GET
    so that the redirect is followed once and every later request -- including the POSTs
    that fetch cDNA, which a redirect would strip of their body -- goes straight to the
    host it lands on.  That host must report exactly ``release``.  The alias is asked
    first, once, so that a pinned run against an archived release does not depend on
    ``rest.ensembl.org`` at all.  Only if that fails is :data:`SERVER` asked which
    release is current: the current release has no working alias and is served by
    :data:`SERVER`, a later one is not served, and for an earlier one the alias is asked
    again under the full retry policy.

    Raises :class:`ReleaseNotServed` for a release later than the current one, for an
    archive that reports a different release, for an archive whose alias now redirects
    to Ensembl's page on archives -- which is how Ensembl retires one -- and for an archive
    that cannot be reached once the retries are spent, which is an outage or a retirement
    that looks like one; a redirect that ends elsewhere off the REST service is one of
    those.
    """
    if release is None:
        return SERVER
    release = int(release)
    net = {"timeout": timeout, "retries": retries, "retry_wait": retry_wait}
    alias = ARCHIVE_ALIAS % release
    try:
        return _archive(alias, release, **dict(net, retries=0))
    except TRANSIENT + (_NotAListing,):
        pass                        # no alias for the current release, or a blip
    current = fetch_release(**net)
    if release == current:
        return SERVER
    if release > current:
        if current >= REST_LAST_RELEASE:
            raise ReleaseNotServed(
                "Ensembl release %d is not on the REST API: %s serves release %d, the last "
                "release Ensembl publishes there; later releases are published only on "
                "the new Ensembl platform, which has no REST API" % (release, SERVER, current))
        raise ReleaseNotServed("Ensembl release %d does not exist yet: %s serves %d, the "
                               "current release" % (release, SERVER, current))
    try:
        return _archive(alias, release, **net)
    except TRANSIENT + (_NotAListing,) as e:
        raise ReleaseNotServed(
            "could not reach the REST archive for Ensembl release %d at %s (%s). It is "
            "either temporarily down or retired; Ensembl keeps REST archives for recent "
            "releases only." % (release, alias, e)) from e


def fetch_cdna_batch(ids, versions=None, **retry):
    """``{id: cdna}`` for Ensembl transcript ids, :data:`MAX_POST_IDS` per request.

    Versions are stripped from the ids.  Results are matched to ids by the ``query``
    field Ensembl echoes back (``id`` if it is absent), never by position.  An id for
    which Ensembl returns nothing is absent from the result; the caller decides whether
    that is an error.  ``versions``, a dict, receives the versioned id of each sequence
    returned, when the answer gives one.  ``retry`` goes to :func:`request`.
    """
    want = list(dict.fromkeys(t.split(".")[0] for t in ids))
    out = {}
    for i in range(0, len(want), MAX_POST_IDS):
        chunk = want[i:i + MAX_POST_IDS]
        for item in request_json("/sequence/id?type=cdna", {"ids": chunk}, **retry):
            tid = (item.get("query") or item["id"]).split(".")[0]
            out[tid] = item["seq"]
            if versions is not None:
                full = str(item.get("id") or "")
                if "." not in full and item.get("version") is not None:
                    full = "%s.%s" % (tid, item["version"])
                if "." in full:
                    versions[tid] = full
    return out
