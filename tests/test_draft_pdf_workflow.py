"""The JOSS draft-PDF workflow runs when the paper changes, and not on every push.

Read with regular expressions, as the other workflow and metadata tests read their files.
"""
import pathlib
import re

WORKFLOW = pathlib.Path(__file__).resolve().parents[1] / ".github" / "workflows" / "draft-pdf.yml"
PATHS = ["paper.md", "paper.bib", ".github/workflows/draft-pdf.yml"]


def _block(text, key):
    """The lines indented under ``key:`` (two spaces in), up to the next key at its level."""
    m = re.search(r"^  %s:\n((?:    .*\n)+)" % re.escape(key), text, re.M)
    assert m, key
    return m.group(1)


def test_it_runs_on_the_paper_alone():
    text = WORKFLOW.read_text()
    push, pr = _block(text, "push"), _block(text, "pull_request")
    assert re.search(r"^    branches: \[main\]$", push, re.M)
    for block in (push, pr):
        assert re.findall(r"^      - (\S+)$", block, re.M) == PATHS
    assert re.search(r"^  workflow_dispatch:$", text, re.M)


def test_it_uses_the_joss_action_and_keeps_the_pdf():
    text = WORKFLOW.read_text()
    assert "uses: actions/checkout@v5" in text
    assert "uses: openjournals/openjournals-draft-action@master" in text
    assert re.search(r"^\s+journal: joss$", text, re.M)
    assert re.search(r"^\s+paper-path: paper\.md$", text, re.M)
    assert "uses: actions/upload-artifact@v4" in text
    assert re.search(r"^\s+path: paper\.pdf$", text, re.M)
