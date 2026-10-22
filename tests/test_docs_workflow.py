"""The documentation site is built strictly on every change to it, and deployed from main.

``.github/workflows/docs.yml`` builds the site with ``mkdocs build --strict`` -- a broken
link or a page missing from the nav fails the build -- and uploads it as the Pages
artifact; a second job deploys it, only from a push to (or a manual run on) ``main`` of the
public repository, so a pull request or a private mirror's push builds and stops there.
Only that job may write to Pages, and only it is in the ``pages`` concurrency group: put on
the whole workflow, the group would hold a pull request's build in the queue of a pending
deployment, and cancel one.  Nothing is pushed to a ``gh-pages`` branch.  Read with regular
expressions, as the other workflow tests read their files.
"""
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "docs.yml"
#: the sdist (the folder with PKG-INFO) leaves .github/ out; in a checkout these tests run,
#: and a workflow that is not there fails them rather than skipping them
pytestmark = pytest.mark.skipif(
    (ROOT / "PKG-INFO").exists() and not WORKFLOW.exists(),
    reason=".github/workflows/docs.yml is not here (the sdist leaves .github/ out)")
PATHS = ["docs/**", "mkdocs.yml", "README.md", "scripts/make_docs_index.py",
         ".github/workflows/docs.yml"]
DEPLOY_IF = ("github.repository == 'charliekim97/isoform-dominance-pipeline' && "
             "github.ref == 'refs/heads/main' && "
             "(github.event_name == 'push' || github.event_name == 'workflow_dispatch')")


def _text():
    return WORKFLOW.read_text(encoding="utf-8")


def _block(text, key, indent):
    """The lines indented under ``key:`` at ``indent`` spaces, up to the next key at that
    level or above."""
    m = re.search(r"^%s%s:[ \t]*\n((?:(?:%s .*)?\n)+)" % (" " * indent, re.escape(key),
                                                          " " * indent), text, re.M)
    assert m, key
    return m.group(1)


def _job(name):
    return _block(_block(_text(), "jobs", 0), name, 2)


def _paths(block):
    return [p.strip("\"'") for p in re.findall(r"^\s+- (\S+)\s*$", block, re.M)]


def test_it_runs_on_the_docs_alone_and_by_hand():
    on = _block(_text(), "on", 0)
    push, pr = _block(on, "push", 2), _block(on, "pull_request", 2)
    assert re.search(r"^    branches: \[main\]$", push, re.M)
    for block in (push, pr):
        assert _paths(_block(block, "paths", 4)) == PATHS
    assert re.search(r"^  workflow_dispatch:", on, re.M)


def test_the_build_is_strict_and_uploads_the_site():
    build = _job("build")
    assert re.search(r"^    permissions:\n      contents: read$", build, re.M)
    assert "uses: actions/checkout@v5" in build
    assert "uses: actions/setup-python@v6" in build
    assert re.search(r"pip install \"mkdocs==1\.6\.1\"$", build, re.M)
    assert re.search(r"run: mkdocs build --strict --site-dir site$", build, re.M)
    upload = build[build.index("uses: actions/upload-pages-artifact@v5"):]
    assert re.search(r"^\s+with:\n\s+path: site$", upload, re.M)


def test_only_a_push_to_main_of_the_public_repository_deploys():
    deploy = _job("deploy")
    assert re.search(r"^    needs: build$", deploy, re.M)
    assert re.search(r"^    if: (.+)$", deploy, re.M).group(1) == DEPLOY_IF
    assert re.search(r"^    permissions:\n      pages: write\n      id-token: write$", deploy,
                     re.M)
    env = _block(deploy, "environment", 4)
    assert re.search(r"^      name: github-pages$", env, re.M)
    assert re.search(r"^      url: \$\{\{ steps\.deployment\.outputs\.page_url \}\}$", env,
                     re.M)
    assert re.search(r"^      - id: deployment\n        uses: actions/deploy-pages@v5$", deploy,
                     re.M)


def test_only_the_deploy_job_is_in_the_pages_concurrency_group():
    text = _text()
    deploy = _job("deploy")
    assert re.search(r"^    concurrency:\n      group: pages\n      cancel-in-progress: false$",
                     deploy, re.M)
    assert len(re.findall(r"^\s*concurrency:", text, re.M)) == 1
    assert "concurrency" not in _job("build")


def test_nothing_is_pushed_to_a_branch():
    text = _text()
    for word in ("gh-deploy", "gh-pages", "git push", "contents: write"):
        assert word not in text, word
    assert len(re.findall(r"^\s+pages: write$", text, re.M)) == 1


def test_the_site_has_the_four_pages_and_is_not_committed():
    nav = re.findall(r"^\s+- .*?:\s*(\S+\.md)\s*$",
                     (ROOT / "mkdocs.yml").read_text(encoding="utf-8"), re.M)
    assert nav == ["index.md", "tutorial_NTRK2.md", "example_genes.md", "api.md"]
    for page in nav:
        assert (ROOT / "docs" / page).is_file(), page
    assert "site/" in (ROOT / ".gitignore").read_text(encoding="utf-8").split()
