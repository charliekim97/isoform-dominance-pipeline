"""Only a redirect to Ensembl's page on archives says a REST archive is retired.

2.4.1 took any redirect that ended outside the REST service for a retired archive and did
not try again, so a run that met a maintenance page reported a retired archive.  On
2026-10-02 (``/info/ping``) the retired archives ended on the same page and the live ones on
their own date-named host:

    e100 -> apr2020.rest.ensembl.org -> www.ensembl.org/help/articles/archives   (HTML)
    e103 -> feb2021.rest.ensembl.org -> www.ensembl.org/help/articles/archives   (HTML)
    e104 -> may2021.rest.ensembl.org -> www.ensembl.org/help/articles/archives   (HTML)
    e105 -> dec2021.rest.ensembl.org   {"ping":1}
    e110 -> jul2023.rest.ensembl.org   {"ping":1}
    e115 -> sep2025.rest.ensembl.org   {"ping":1}

That page, on ``www.ensembl.org`` or ``ensembl.org``, is retirement.  Any other end outside
the REST service is retried like an outage, and named as one when the retries run out.
"""
import time
import urllib.request

import pytest

from isoform_dominance import ensembl
from test_ensembl_release import ALIAS, HOST, OLD, Ensembl, _Resp


class _Elsewhere(Ensembl):
    """The release-``OLD`` alias ends at ``page`` for its first ``times`` requests."""

    def __init__(self, page, times, **kw):
        super().__init__(**kw)
        self.page, self.left = page, times

    def __call__(self, req, timeout=None, **kw):
        if req.full_url.startswith(ALIAS % OLD) and self.left:
            self.left -= 1
            self.calls.append((req.get_method(), ALIAS % OLD, "/info/data"))
            return _Resp("<html>somewhere else</html>", self.page)
        return super().__call__(req, timeout=timeout, **kw)


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
def elsewhere(monkeypatch, sleeps):
    def _serve(page, times):
        fake = _Elsewhere(page, times)
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        return fake
    return _serve


@pytest.mark.parametrize("page", [
    "https://www.ensembl.org/help/articles/archives",
    "https://ensembl.org/help/articles/archives",
    "http://www.ensembl.org/help/articles/archives?from=e104",
    "https://www.ensembl.org/help/articles/archives/retired",
])
def test_the_archives_page_is_retirement_at_once(elsewhere, sleeps, page):
    fake = elsewhere(page, 99)
    with pytest.raises(ensembl.ReleaseNotServed, match="release 110 is retired"):
        ensembl.resolve_server(OLD)
    assert fake.hosts() == [ALIAS % OLD] and sleeps == []


@pytest.mark.parametrize("page", [
    "https://www.ensembl.org/info/about/maintenance.html",     # not the archives page
    "https://status.ensembl.org/",                              # another host
    "https://www.ensembl.org.example.com/help/articles/archives",
    "https://www.ensembl.org/help/articles/archive-of-old-news",
])
def test_any_other_end_off_the_service_is_tried_again(elsewhere, sleeps, page):
    fake = elsewhere(page, 2)                   # then the archive answers
    assert ensembl.resolve_server(OLD, retries=3, retry_wait=1.0) == HOST[OLD]
    # one try without retries, then the full policy: the alias ends there once more, and
    # after one wait it answers
    assert sleeps == [1.0]
    assert fake.hosts().count(ALIAS % OLD) == 3


def test_an_end_off_the_service_that_stays_is_down_or_retired(elsewhere, sleeps):
    elsewhere("https://www.ensembl.org/info/about/maintenance.html", 99)
    with pytest.raises(ensembl.ReleaseNotServed, match="temporarily down or retired") as e:
        ensembl.resolve_server(OLD, retries=2, retry_wait=1.0)
    assert "maintenance.html" in str(e.value)
    assert sleeps == [1.0, 2.0]


def test_the_date_named_host_is_the_archive(serve, sleeps):
    fake = serve()
    assert ensembl.resolve_server(OLD) == HOST[OLD]
    assert fake.hosts() == [ALIAS % OLD] and sleeps == []
