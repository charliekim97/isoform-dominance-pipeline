"""``paper.md`` against what JOSS asks of a paper, and what this one must keep.

JOSS asks for 750 to 1,750 words and for these sections, in this order: Summary, Statement of
need, State of the field, Software design, Research impact statement, AI usage disclosure,
Acknowledgements, References.  The AI usage disclosure is an HTML comment until it is restored
before submission, at about 150 words, so while it is one the paper stays at or under 1,600
words, and with the section restored at or under 1,750.

The count: the YAML header, HTML comments and heading lines taken out, ``*`` and ``_`` taken
out, everything from ``# References`` on dropped, and the rest split on white space.
"""
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper.md"
BIB = ROOT / "paper.bib"

#: the headings, in order; the AI usage disclosure goes before the acknowledgements once it
#: is restored
HEADINGS = ["Summary", "Statement of need", "State of the field", "Software design",
            "Research impact statement", "Acknowledgements", "References"]
AI_USAGE = "AI usage disclosure"
#: the most words with the AI usage disclosure an HTML comment, and with it restored
LIMIT_WITHOUT_AI_USAGE, LIMIT = 1600, 1750
#: the fewest words JOSS asks for
FLOOR = 750
#: the 23 references the paper cites
KEYS = {"anders2012", "ferrerbonsoms2022", "frankish2021", "gilis2022", "gse137619",
        "gse228458", "harris2020", "hiller2009", "hunter2007", "liptak1958", "nowicka2016",
        "patro2017", "sarkar2020", "soneson2016", "stouffer1949", "vanelteren1960",
        "virtanen2020", "vittingseerup2019", "wilcoxon1945", "yang2026", "yates2015",
        "zheng2022", "zhu2019"}
#: the reference the HTML comment in "Research impact statement" says may go, with the
#: paragraph that cites it, if that study has not been accepted by submission
REMOVABLE = ("yang2026", "*Use in research.*")


def text_of(path=PAPER):
    return path.read_text(encoding="utf-8")


def without_header(text):
    return re.sub(r"\A---\n.*?\n---\n", "", text, flags=re.S)


def body(text):
    """The text the words are counted in, before its heading lines are taken out."""
    text = without_header(text)
    text = text[:text.index("\n# References")]
    return re.sub(r"<!--.*?-->", "", text, flags=re.S)


def word_count(text):
    lines = [line for line in body(text).split("\n") if not line.startswith("#")]
    return len(re.sub(r"[*_]", "", "\n".join(lines)).split())


def headings(text):
    return re.findall(r"^# (.+?)\s*$", without_header(text), re.M)


def cited(text):
    """The citation keys of the text a reader sees: outside HTML comments, before the
    references."""
    return set(re.findall(r"(?<![\w.])@([A-Za-z][\w-]*)", body(text)))


def ai_usage_restored(text):
    return AI_USAGE in headings(text)


def test_the_count_is_the_one_described():
    text = ("---\ntitle: 'x'\n---\n\n# Summary\n\nOne *two* three_four.\n<!-- five six -->\n\n"
            "# Statement of need\n\nseven `eight`\n\n# References\n\nnine ten\n")
    assert word_count(text) == 5


def test_the_paper_is_as_long_as_joss_asks():
    text = text_of()
    n = word_count(text)
    assert n >= FLOOR
    if ai_usage_restored(text):
        assert n <= LIMIT, n
    else:
        # the disclosure's place is held by an HTML comment, which the count leaves out
        assert any(AI_USAGE in c for c in re.findall(r"<!--.*?-->", text, flags=re.S))
        assert n <= LIMIT_WITHOUT_AI_USAGE, n


def test_the_sections_are_joss_sections_in_its_order():
    got = headings(text_of())
    want = list(HEADINGS)
    if AI_USAGE in got:
        want.insert(want.index("Acknowledgements"), AI_USAGE)
    assert got == want


def test_the_paper_cites_its_references():
    text = text_of()
    keys = cited(text)
    key, paragraph = REMOVABLE
    if keys == KEYS - {key}:
        assert paragraph not in text
    else:
        assert keys == KEYS
    bib = set(re.findall(r"^@\w+\{([^,\s]+),", text_of(BIB), re.M))
    assert keys <= bib


@pytest.mark.parametrize("change", ["a word too many", "a key gone", "a heading renamed",
                                    "a section moved"])
def test_each_check_can_fail(change):
    text = text_of()
    if change == "a word too many":
        n = word_count(text)
        limit = LIMIT if ai_usage_restored(text) else LIMIT_WITHOUT_AI_USAGE
        filler = " ".join(["word"] * (limit + 1 - n))
        text = text.replace("\n# Acknowledgements", "\n%s\n\n# Acknowledgements" % filler, 1)
        assert word_count(text) == limit + 1
    elif change == "a key gone":
        text = text.replace("[@hunter2007]", "")
        assert cited(text) not in (KEYS, KEYS - {REMOVABLE[0]})
    elif change == "a heading renamed":
        text = text.replace("\n# Research impact statement\n", "\n# Research impact\n")
        assert headings(text) != HEADINGS
    else:
        sd = text.index("\n# Software design\n")
        sotf = text.index("\n# State of the field\n")
        text = text[:sotf] + text[sd:text.index("\n# Research impact")] + text[sotf:sd] + \
            text[text.index("\n# Research impact"):]
        assert headings(text)[2:4] == ["Software design", "State of the field"]
