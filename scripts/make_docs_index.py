"""Write docs/index.md, the documentation site's front page, from README.md.

The README is PyPI's project page too, where a relative link resolves against pypi.org and
breaks, so it names the figure and the other documentation pages by URL: the figure by its
raw file on the main branch, a page by its view on GitHub.  The site holds those files beside
index.md, so here each such URL becomes the file's name again, and nothing else changes.
Run it after editing the README; ``tests/test_docs_index.py`` fails until the two agree::

    python scripts/make_docs_index.py
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
REPO = "charliekim97/isoform-dominance-pipeline"
#: a file of docs/ as the README links it: an image by its raw file, a page by its view
RAW = "https://raw.githubusercontent.com/%s/main/docs/" % REPO
BLOB = "https://github.com/%s/blob/main/docs/" % REPO


def index_text(readme):
    """The README's text with each link to a file of docs/ made relative to docs/: a
    ``RAW`` URL of any file there, and a ``BLOB`` URL of a Markdown page there."""
    text = re.sub(r"\]\(%s([\w.-]+)\)" % re.escape(RAW), r"](\1)", readme)
    return re.sub(r"\]\(%s([\w.-]+\.md)\)" % re.escape(BLOB), r"](\1)", text)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    out = pathlib.Path(argv[0]) if argv else ROOT / "docs" / "index.md"
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write(index_text(readme))
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
