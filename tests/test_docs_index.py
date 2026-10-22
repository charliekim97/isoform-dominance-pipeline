"""The README links by absolute URL, and the docs' front page is the README.

PyPI shows the README as the project page, where a relative link resolves against
pypi.org and breaks: the figure and the three documentation pages were four such links.
The README now names them by URL, and ``scripts/make_docs_index.py`` writes
``docs/index.md`` from it, those URLs made relative again for the site, which holds the
files beside the page.  The committed page must be what the script writes now, so a README
edited without regenerating it fails here.  Text, not bytes, is compared: a Windows
checkout may turn the page's line ends into CRLF.
"""
import importlib.util
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
INDEX = ROOT / "docs" / "index.md"


def _script():
    spec = importlib.util.spec_from_file_location(
        "_make_docs_index", ROOT / "scripts" / "make_docs_index.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _prose(text):
    """``text`` without its fenced code blocks and code spans: what is rendered as links."""
    text = re.sub(r"^```.*?^```[ \t]*$", "", text, flags=re.M | re.S)
    return re.sub(r"`[^`\n]*`", "", text)


def link_targets(text):
    """Every target of an inline Markdown link or image, a reference definition, and an
    HTML ``src`` or ``href``, outside code."""
    prose = _prose(text)
    return (re.findall(r"\]\(\s*<?([^)\s>]*)", prose)
            + re.findall(r"^\s{0,3}\[[^\]]+\]:\s*<?(\S+?)>?(?:\s|$)", prose, re.M)
            + re.findall(r"\b(?:src|href)\s*=\s*[\"']([^\"']*)", prose))


def relative(targets):
    return [t for t in targets if not re.match(r"(?:https?:|mailto:|#)", t)]


def test_the_readme_has_no_relative_link():
    targets = link_targets(README.read_text(encoding="utf-8"))
    assert len(targets) > 10, targets
    assert relative(targets) == []


def test_the_check_finds_a_relative_link():
    # the check above, on the README as it was in 2.8.0, finds its four relative links
    old = ("![example output](docs/example_output.png)\n[Statistical notes](#statistical-notes)\n"
           "[NTRK2 walkthrough](docs/tutorial_NTRK2.md)\n[ref]: docs/api.md\n"
           "<img src=\"docs/example_output.svg\">\n```\n[not a link](inside/code.md)\n```\n")
    assert relative(link_targets(old)) == ["docs/example_output.png", "docs/tutorial_NTRK2.md",
                                           "docs/api.md", "docs/example_output.svg"]


def test_the_front_page_is_the_readme_as_the_script_writes_it():
    made = _script().index_text(README.read_text(encoding="utf-8"))
    assert INDEX.read_text(encoding="utf-8") == made


def test_the_front_page_links_the_files_beside_it():
    text = INDEX.read_text(encoding="utf-8")
    targets = relative(link_targets(text))
    assert sorted(targets) == ["api.md", "example_genes.md", "example_output.png",
                               "tutorial_NTRK2.md"]
    for t in targets:
        assert (INDEX.parent / t).is_file(), t


def test_the_script_changes_only_the_two_kinds_of_url():
    script = _script()
    readme = README.read_text(encoding="utf-8")
    made = script.index_text(readme)
    back = made
    for page in ("tutorial_NTRK2.md", "example_genes.md", "api.md"):
        back = back.replace("](%s)" % page, "](%s%s)" % (script.BLOB, page))
    back = back.replace("](example_output.png)", "](%sexample_output.png)" % script.RAW)
    assert back == readme
    # another file of the repository keeps its URL: the site does not hold it
    other = "[x](%sexample_output_stats.csv) [y](%s../README.md)" % (script.BLOB, script.BLOB)
    assert script.index_text(other) == other
