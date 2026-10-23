"""The NTRK2 walkthrough and the example-genes page show what the commands print, offline.

Both pages run on the inputs saved in ``example/saved/``, from the repository root, with
every socket refused.  A page has two kinds of block:

- ```` ```bash ````: one ``isoform-dominance`` command (``\\`` continues it), then its
  standard output, each line after ``# ``, and ``# exit status N``; run here with
  ``cli.main`` and compared line by line, with the exit status;
- ```` ```python ````: a few lines that use the package, then their standard output, each
  line after ``# ``; run here with ``exec`` and compared line by line.

``annotate`` needs a GTF or Ensembl's REST API.  An ``annotate --gene G --gtf
gencode.vN.annotation.gtf.gz`` block runs with ``annotate.run`` returning
``example/saved/G.gencodeN.config.json``, the config that command wrote, so what it prints
is what the command printed; the GTF itself is not in the repository.  A quotation on a
page (``> ``) is what a command on that page wrote on standard error, word for word.  The
example-genes table is recomputed from the saved inputs.
"""
import copy
import json
import pathlib
import re
import shlex
import socket

import pytest

from isoform_dominance import annotate, cli

ROOT = pathlib.Path(__file__).resolve().parents[1]
SAVED = ROOT / "example" / "saved"
#: the pages, and the most lines each may run to
PAGES = {"docs/tutorial_NTRK2.md": 128, "docs/example_genes.md": 92}


def text_of(page):
    return (ROOT / page).read_text(encoding="utf-8")


def fenced(text):
    """``[(language, body)]`` of the fenced code blocks of a Markdown text."""
    return re.findall(r"^```([\w-]*)[ \t]*\n(.*?)^```[ \t]*$", text, re.M | re.S)


def bash_block(body):
    """``(argv, the standard output shown, the exit status shown)`` of a bash block.  Every
    line of the command but its last ends with `` \\``, as a shell needs to read one command."""
    lines = body.rstrip("\n").split("\n")
    first = next(i for i, x in enumerate(lines) if x.startswith("#"))
    assert first, body
    assert all(x.endswith(" \\") for x in lines[:first - 1]), lines[:first]
    assert not lines[first - 1].endswith("\\"), lines[:first]
    command = " ".join(x[:-1].strip() if x.endswith("\\") else x.strip()
                       for x in lines[:first])
    shown = lines[first:]
    assert all(x == "#" or x.startswith("# ") for x in shown), shown
    m = re.fullmatch(r"# exit status (\d+)", shown[-1])
    assert m, shown[-1]
    argv = shlex.split(command)
    assert argv[0] == "isoform-dominance", command
    return argv[1:], [x[2:] for x in shown[:-1]], int(m.group(1))


def python_block(body):
    """``(code, the standard output shown)`` of a python block: the output is the run of
    ``# `` lines that ends it."""
    lines = body.rstrip("\n").split("\n")
    start = len(lines)
    while start and (lines[start - 1] == "#" or lines[start - 1].startswith("# ")):
        start -= 1
    assert 0 < start < len(lines), body
    return "\n".join(lines[:start]) + "\n", [x[2:] for x in lines[start:]]


@pytest.fixture
def offline(monkeypatch):
    """No connection can be made, and the working directory is the repository's root."""
    def refuse(*a, **k):
        raise AssertionError("a documented example tried to reach the network")
    monkeypatch.setattr(socket, "socket", refuse)
    monkeypatch.setattr(socket, "create_connection", refuse)
    monkeypatch.setattr(socket, "getaddrinfo", refuse)
    monkeypatch.chdir(ROOT)


def run_command(argv, monkeypatch, capsys):
    """``(stdout lines, stderr, exit status)`` of a documented command, run in-process."""
    if argv[0] == "annotate":
        gene = argv[argv.index("--gene") + 1]
        gtf = argv[argv.index("--gtf") + 1] if "--gtf" in argv else ""
        m = re.fullmatch(r"gencode\.v(\d+)\.annotation\.gtf\.gz", gtf)
        assert m, "an annotate example names its GENCODE GTF: %s" % argv
        with open(SAVED / ("%s.gencode%s.config.json" % (gene, m.group(1))),
                  encoding="utf-8") as f:
            wrote = json.load(f)

        def run(g, out, **kw):
            assert (g, kw["gtf"]) == (gene, gtf)
            return copy.deepcopy(wrote)
        monkeypatch.setattr(annotate, "run", run)
    capsys.readouterr()
    code = cli.main(argv)
    out, err = capsys.readouterr()
    return out.rstrip("\n").split("\n") if out else [], err, code


def run_python(code, capsys):
    capsys.readouterr()
    exec(compile(code, "<docs>", "exec"), {"__name__": "__docs__"})
    out = capsys.readouterr().out
    return out.rstrip("\n").split("\n") if out else []


BLOCKS = [(page, i, lang, body) for page in PAGES
          for i, (lang, body) in enumerate(fenced(text_of(page)))]


def test_every_block_is_one_of_the_two_kinds():
    assert {lang for _, _, lang, _ in BLOCKS} == {"bash", "python"}
    for _, _, lang, body in BLOCKS:
        bash_block(body) if lang == "bash" else python_block(body)


@pytest.mark.parametrize("page", sorted(PAGES))
def test_every_code_block_is_one_the_tests_run(page):
    # a block the tests do not find is a block nothing checks: a ~~~ fence, a fence with
    # more to its info string than a language, or an indented block
    text = text_of(page)
    lines = text.split("\n")
    assert not [x for x in lines if x.lstrip().startswith("~~~")]
    assert len([x for x in lines if x.startswith("```")]) == 2 * len(fenced(text))
    outside = re.sub(r"^```.*?^```[ \t]*$", "", text, flags=re.M | re.S).split("\n")
    assert [x for x in outside if x.startswith(("    ", "\t")) and x.strip()] == []


@pytest.mark.parametrize("page,i,lang,body", BLOCKS,
                         ids=["%s#%d" % (p.split("/")[-1], i) for p, i, _, _ in BLOCKS])
def test_each_block_shows_what_it_prints(offline, monkeypatch, capsys, page, i, lang, body):
    if lang == "bash":
        argv, shown, status = bash_block(body)
        out, _, code = run_command(argv, monkeypatch, capsys)
        assert (out, code) == (shown, status)
    else:
        code, shown = python_block(body)
        assert run_python(code, capsys) == shown


def _quotes(text):
    """Each quotation of a page, its lines joined and its white space made single."""
    quotes, current = [], []
    for line in text.split("\n") + [""]:
        if line.startswith(">"):
            current.append(line[1:].strip())
        elif current:
            quotes.append(" ".join(" ".join(current).split()))
            current = []
    return quotes


@pytest.mark.parametrize("page", sorted(PAGES))
def test_each_quotation_is_what_a_command_on_the_page_wrote(offline, monkeypatch, capsys,
                                                            page):
    quotes = _quotes(text_of(page))
    assert quotes
    said = []
    for _, _, _, body in [b for b in BLOCKS if b[0] == page and b[2] == "bash"]:
        said.append(run_command(bash_block(body)[0], monkeypatch, capsys)[1])
    said = " ".join(" ".join(said).split())
    for q in quotes:
        assert q in said, q


@pytest.mark.parametrize("page", sorted(PAGES))
def test_each_page_says_where_it_runs(page):
    text = text_of(page)
    head = text[:text.index("\n## ")] if "\n## " in text else text
    assert "clone" in head and "`example/saved/`" in head
    assert len(text.rstrip("\n").split("\n")) <= PAGES[page]


# --------------------------------------------------------------------------- #
# the example-genes table
# --------------------------------------------------------------------------- #
GENES = ("LEPR", "NTRK2", "NTRK3", "FLT1")
RELEASES = {44: 110, 50: 116}


def table(text):
    """The rows of the page's table, as lists of cells."""
    rows = [line for line in text.split("\n") if line.startswith("| ")]
    cells = [[c.strip() for c in re.split(r"(?<!\\)\|", r)[1:-1]] for r in rows]
    return cells[0], [r for r in cells[1:] if not set("".join(r)) <= set("-:")]


def _cell(x, fmt):
    return "n/a" if x is None else fmt % x


def row(gene, n, monkeypatch, capsys):
    """The row the saved inputs of ``gene`` at GENCODE ``n`` give."""
    base = "example/saved/%s.gencode%d" % (gene, n)
    argv = ["identifiability", "--config", base + ".config.json",
            "--inputs", base + ".inputs.json"]
    out, _, code = run_command(argv + ["--json"], monkeypatch, capsys)
    assert code == 0
    res = json.loads("\n".join(out))
    with open(ROOT / (base + ".config.json"), encoding="utf-8") as f:
        groups = json.load(f)["groups"]
    a, b = res["primary_comparison"]
    c = res["contrast"]
    _, _, status = run_command(argv + ["--min-log2fc", "0.5"], monkeypatch, capsys)
    return [gene, "%d (%d)" % (n, RELEASES[n]),
            "`%s` (%d) vs `%s` (%d)" % (a, len(groups[a]), b, len(groups[b])),
            "`%s`" % res["verdict"], _cell(c["min_resolvable_log2fc"], "%.3f"),
            # the CLI prints an infinite conditioning factor as inf
            "%.2f" % (c["conditioning_factor"] if c["conditioning_factor"] is not None
                      else float("inf")), str(status)]


HEAD = ["Gene", "GENCODE (Ensembl)", "Proposed comparison (transcripts)", "Verdict",
        "Contrast min \\|log2FC\\|", "Contrast conditioning",
        "Exit status, `--min-log2fc 0.5`"]


def test_the_table_is_what_the_saved_inputs_give(offline, monkeypatch, capsys):
    head, rows = table(text_of("docs/example_genes.md"))
    assert head == HEAD                     # the columns are in the order row() gives them
    assert rows == [row(g, n, monkeypatch, capsys) for g in GENES for n in (44, 50)]
