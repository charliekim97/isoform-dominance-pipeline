"""The weekly job that holds the GENCODE 50 files to REST 116, run against an offline fake.

The step script is taken from the workflow file itself (``ensembl-nightly.yml``, job
``gtf-parity``) and run with an ``isoform-dominance`` on PATH whose Ensembl is
:class:`LiveRest`: REST at release 116 for the genes of ``tests/data/gencode_mini``, from
the recorded lookups and the transcript FASTA.  An outage is skipped, a retired archive
warned of, and an answer that moved warned of, never failed.
"""
import email.message
import gzip
import io
import json
import os
import pathlib
import subprocess
import sys
from urllib.error import HTTPError

import pytest

from isoform_dominance import ensembl
from test_nightly_workflow import step_script

TESTS = os.path.dirname(os.path.abspath(__file__))
DATA = pathlib.Path(TESTS) / "data" / "gencode_mini"
GTF = str(DATA / "gencode.v50.mini.gtf.gz")
FASTA = str(DATA / "gencode.v50.mini.transcripts.fa.gz")
REST = json.loads(gzip.decompress((DATA / "rest116_mini.json.gz").read_bytes()))
RETIRED_PAGE = "https://www.ensembl.org/help/articles/archives"
RENAMED = "ENST00009999999"

pytestmark = pytest.mark.skipif(os.name == "nt", reason="runs a POSIX shim on PATH")


def _fasta():
    out, cur = {}, None
    with gzip.open(FASTA, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split("|")[0]
                out[cur.split(".")[0]] = [cur, []]
            else:
                out[cur.split(".")[0]][1].append(line.strip())
    return {t: (v, "".join(s)) for t, (v, s) in out.items()}


class _Resp(io.BytesIO):
    status = 200

    def __init__(self, payload, url):
        super().__init__(payload.encode() if isinstance(payload, str) else payload)
        self.headers = email.message.Message()
        self.url = url

    def geturl(self):
        return self.url


class LiveRest:
    """rest.ensembl.org serving release 116 for the extract's genes; ``e116``'s alias
    answers 503, as it did live.  ``down``: every host answers 503.  ``retired``: the
    alias redirects to Ensembl's page on archives.  ``moved``: LEPR's lookup has lost a
    protein-coding transcript.  ``renamed``: LepRb has another id, which no GTF holds.
    ``unversioned``: cDNA comes without its version.  ``cdna_down``: cDNA answers 503."""

    def __init__(self, down=False, retired=False, moved=False, unversioned=False,
                 renamed=False, cdna_down=False):
        self.down, self.retired, self.moved = down, retired, moved
        self.unversioned, self.renamed, self.cdna_down = unversioned, renamed, cdna_down
        self.seqs = _fasta()
        if renamed:
            self.seqs[RENAMED] = ("%s.1" % RENAMED, self.seqs["ENST00000349533"][1])

    def _lookup(self, gid):
        g = REST["lookups"][gid]
        if self.moved and g["display_name"] == "LEPR":
            drop = next(t["id"] for t in g["Transcript"] if t["biotype"] == "protein_coding"
                        and t["id"] != "ENST00000349533")
            g = dict(g, Transcript=[t for t in g["Transcript"] if t["id"] != drop])
        if self.renamed and g["display_name"] == "LEPR":
            g = dict(g, canonical_transcript="%s.1" % RENAMED, Transcript=[
                dict(t, id=RENAMED, version=1) if t["id"] == "ENST00000349533" else t
                for t in g["Transcript"]])
        return g

    def __call__(self, req, timeout=None, **_):
        url, method = req.full_url, req.get_method()
        host = "/".join(url.split("/")[:3])
        path = url[len(host):]
        if self.down:
            raise HTTPError(url, 503, "status 503", email.message.Message(), io.BytesIO())
        if host == ensembl.ARCHIVE_ALIAS % 116:
            if self.retired:
                return _Resp("<html>Ensembl archives</html>", RETIRED_PAGE)
            raise HTTPError(url, 503, "status 503", email.message.Message(), io.BytesIO())
        assert host == ensembl.SERVER, url
        parts = path.split("?")[0].split("/")
        if path.startswith("/info/data"):
            return _Resp(json.dumps({"releases": [116]}), url)
        names = {}
        for g in REST["lookups"].values():
            names.setdefault(g["display_name"], []).append(g["id"])
        if path.startswith("/lookup/symbol/"):
            gid = REST["lookup_symbol"].get(parts[4]) or names[parts[4]][0]
            return _Resp(json.dumps(self._lookup(gid)), url)
        if path.startswith("/xrefs/symbol/"):
            return _Resp(json.dumps([{"type": "gene", "id": i}
                                     for i in names.get(parts[4], [])]), url)
        if path.startswith("/lookup/id/"):
            return _Resp(json.dumps(self._lookup(parts[3])), url)
        body = json.loads(req.data) if req.data else None
        if method == "POST" and path.startswith("/lookup/id"):
            return _Resp(json.dumps({i: self._lookup(i) for i in body["ids"]
                                     if i in REST["lookups"]}), url)
        if method == "POST" and path.startswith("/sequence/id"):
            if self.cdna_down:
                raise HTTPError(url, 503, "status 503", email.message.Message(), io.BytesIO())
            return _Resp(json.dumps([
                dict({"query": i, "seq": self.seqs[i][1]},
                     **({"id": i} if self.unversioned else
                        {"id": self.seqs[i][0], "version": int(self.seqs[i][0].split(".")[1])}))
                for i in body["ids"] if i in self.seqs]), url)
        raise AssertionError("unexpected request %s %s" % (method, url))


def _env(tmp_path, **fake):
    """The environment of the step, with an ``isoform-dominance`` on PATH whose Ensembl is
    ``LiveRest(**fake)``."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    shim = bindir / "isoform-dominance"
    shim.write_text(
        "#!%s\n"
        "import sys, time, urllib.request\n"
        "sys.path.insert(0, %r)\n"
        "from test_gtf_live_workflow import LiveRest\n"
        "time.sleep = lambda s: None\n"
        "urllib.request.urlopen = LiveRest(**%r)\n"
        "from isoform_dominance import cli\n"
        "sys.exit(cli.main(sys.argv[1:]))\n" % (sys.executable, TESTS, fake))
    shim.chmod(0o755)
    return dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ["PATH"], GTF=GTF,
                FASTA=FASTA)


def _run(tmp_path, env=None, **fake):
    script = step_script("parity").replace("/tmp/", str(tmp_path) + "/")
    return subprocess.run([sys.executable, "-"], input=script,
                          env=env or _env(tmp_path, **fake), capture_output=True, text=True)


def test_rest_116_and_the_files_agree(tmp_path):
    r = _run(tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "ok: REST 116 and the GENCODE 50 files agree" in r.stdout
    assert "::warning::" not in r.stdout and "::notice::" not in r.stdout


def test_an_outage_is_skipped(tmp_path):
    r = _run(tmp_path, down=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "::notice::annotate could not reach Ensembl" in r.stdout
    assert "::warning::" not in r.stdout


def test_a_retired_archive_is_a_warning(tmp_path):
    r = _run(tmp_path, retired=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "::warning::the REST archive for release 116 is retired" in r.stdout


def test_a_moved_answer_is_a_warning_naming_what_moved(tmp_path):
    """LEPR's lookup lost a transcript, and REST's cDNA lost its version: each is named."""
    r = _run(tmp_path, moved=True, unversioned=True)
    assert r.returncode == 0, r.stdout + r.stderr
    warning = next(x for x in r.stdout.splitlines() if x.startswith("::warning::"))
    assert warning.startswith("::warning::REST 116 and the GENCODE 50 files no longer agree: "
                              "annotate LEPR (")
    assert "identifiability LEPR" in warning and "CD99" not in warning
    assert "REST gave no version for ENST" in warning


def test_a_failure_of_the_files_side_fails(tmp_path):
    """The files make no request: their failing is this package's, not Ensembl's."""
    env = dict(_env(tmp_path), GTF=str(tmp_path / "gone.gtf.gz"))
    r = _run(tmp_path, env=env)
    assert r.returncode == 1 and "gone.gtf.gz" in r.stdout + r.stderr


def test_an_id_no_gtf_holds_is_a_moved_answer_not_a_failure(tmp_path):
    """The files' side runs on the config the GTF gives, never on REST's: an id REST
    names that the GTF does not is a moved answer."""
    r = _run(tmp_path, renamed=True)
    assert r.returncode == 0, r.stdout + r.stderr
    warning = next(x for x in r.stdout.splitlines() if x.startswith("::warning::"))
    assert "annotate LEPR (" in warning and "identifiability LEPR" in warning


def test_an_outage_after_a_moved_answer_still_says_what_moved(tmp_path):
    r = _run(tmp_path, moved=True, cdna_down=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "::warning::REST 116 and the GENCODE 50 files no longer agree: annotate LEPR (" \
        in r.stdout
    assert "::notice::identifiability could not reach Ensembl" in r.stdout
