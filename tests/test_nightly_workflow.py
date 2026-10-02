"""The weekly live check sorts a retired REST archive from an outage.

Its scripts skipped any failure whose output held one of the strings in ``NET``, and one of
them, "Ensembl release not available", begins the message for a retired archive too ("...:
the REST archive for Ensembl release 110 is retired"), so the release-110 target would have
skipped quietly every week once Ensembl retired that archive.  The step scripts below are
taken from the workflow file itself and run against the offline fake of the archive.
"""
import json
import os
import re
import subprocess
import sys
import textwrap

import pytest

WORKFLOW = os.path.join(os.path.dirname(__file__), os.pardir, ".github", "workflows",
                        "ensembl-nightly.yml")
TESTS = os.path.dirname(os.path.abspath(__file__))

pytestmark = pytest.mark.skipif(os.name == "nt", reason="runs a POSIX shim on PATH")


def step_script(name):
    """The Python heredoc of the workflow step called ``name``."""
    lines = open(WORKFLOW).read().splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.strip() == "- name: %s" % name)
    run = next(i for i in range(start, len(lines)) if lines[i].strip() == "run: |")
    body = []
    for ln in lines[run + 1:]:
        if ln.strip() and len(ln) - len(ln.lstrip()) <= len(lines[run]) - len(
                lines[run].lstrip()):
            break
        body.append(ln)
    text = textwrap.dedent("\n".join(body))
    return re.search(r"python - <<'PY'\n(.*?)\nPY", text, re.S).group(1)


def _shim(tmp_path, **fake):
    """An ``isoform-dominance`` on PATH that talks to the offline fake of Ensembl."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    shim = bindir / "isoform-dominance"
    shim.write_text(
        "#!%s\n"
        "import sys, time, urllib.request\n"
        "sys.path.insert(0, %r)\n"
        "from test_ensembl_release import Ensembl\n"
        "time.sleep = lambda s: None\n"
        "urllib.request.urlopen = Ensembl(**%r)\n"
        "from isoform_dominance import cli\n"
        "sys.exit(cli.main(sys.argv[1:]))\n" % (sys.executable, TESTS, fake))
    shim.chmod(0o755)
    return dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ["PATH"],
                RELEASE="110")


def _run(name, env, tmp_path):
    script = step_script(name).replace("/tmp/", str(tmp_path) + "/")
    return subprocess.run([sys.executable, "-"], input=script, env=env,
                          capture_output=True, text=True)


@pytest.mark.parametrize("step", ["annotate", "identifiability"])
def test_a_retired_archive_is_a_warning_not_a_skipped_outage(tmp_path, step):
    env = _shim(tmp_path, alias_retired=True)
    if step == "identifiability":
        (tmp_path / "lepr.json").write_text(json.dumps(
            {"gene": "FAKE", "groups": {"A": ["ENST00000000001"], "B": ["ENST00000000101"]},
             "primary_comparison": ["A", "B"], "ensembl_release": 110}))
    r = _run(step, env, tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "::warning::" in r.stdout and "retired" in r.stdout
    assert "could not reach Ensembl" not in r.stdout


def test_an_archive_that_does_not_answer_is_still_skipped(tmp_path):
    r = _run("annotate", _shim(tmp_path, alias_down=True, server_down=True), tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "::notice::annotate could not reach Ensembl" in r.stdout
    assert "::warning::" not in r.stdout
