"""One title and one line of introduction, wherever the package names itself.

The title is the paper's, and CITATION.cff carries it, because Zenodo takes its record's
title from there.  The line of introduction is pyproject's ``description`` (PyPI's summary),
said the same in the CLI's help, the package docstring, the README and the docs' front page;
Markdown emphasis is taken off before comparing.  The files that are not Python are read
with regular expressions, as tests/test_dependencies.py reads pyproject.toml.
"""
import pathlib
import re

import isoform_dominance
from isoform_dominance import cli

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _one(pattern, path):
    found = re.findall(pattern, (ROOT / path).read_text(encoding="utf-8"), re.M)
    assert len(found) == 1, (path, found)
    return found[0]


def _plain(text):
    """Markdown emphasis off and white space made single."""
    return " ".join(re.sub(r"\*+|(?<!\w)_+|_+(?!\w)", "", text).split())


def _bold_line(path):
    """The first line of a Markdown file that is bold from end to end."""
    return _one(r"\A(?:.*\n)*?\*\*(.+)\*\*[ \t]*$", path)


def test_the_citation_title_is_the_papers():
    paper = _one(r"^title:\s*'(.+)'\s*$", "paper.md")
    assert paper == ("isoform-dominance: deciding whether a functional isoform-class "
                     "comparison is measurable before running it")
    assert _one(r'^title:\s*"(.+)"\s*$', "CITATION.cff") == paper


def test_one_line_of_introduction_everywhere():
    line = _one(r'^description\s*=\s*"(.+)"\s*$', "pyproject.toml")
    assert line == ("Is this functional isoform-class comparison measurable by short reads? "
                    "Group-level identifiability, quantification and testing from bulk RNA-seq")
    first_paragraph = isoform_dominance.__doc__.strip().split("\n\n")[0]
    for where, text in (("cli", cli.build_parser().description),
                        ("__init__", first_paragraph),
                        ("README.md", _bold_line("README.md")),
                        ("docs/index.md", _bold_line("docs/index.md"))):
        assert _plain(text) == line, where
