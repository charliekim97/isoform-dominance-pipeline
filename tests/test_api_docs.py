"""docs/api.md names every public function and class, with the signature the code has.

A public name is a function or class defined in the module itself (``obj.__module__`` is
the module's) whose name has no leading underscore; a name a module imports from another
is the other module's.  Each has one entry line: ``**`module.name(signature)`**`` for a
function, its signature exactly ``str(inspect.signature(f))``, and
``**`module.Name`** (a subclass of `Base`)`` for a class.  None of these functions has an
annotation or a default whose repr varies, so the strings are the same on Python 3.10 to
3.14.  And every ``module.name`` the page writes in code is a name the module has: one
inside a host name (``rest.ensembl.org``) or a URL is no reference.
"""
import importlib
import inspect
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
API = ROOT / "docs" / "api.md"
MODULES = ("annotate", "annotation_files", "identifiability", "extract", "stats",
           "contamination", "io", "index_scope", "ensembl")
#: what a reference in code may also name: the CLI module, whose functions are internal
REFERABLE = MODULES + ("cli",)
ENTRY = re.compile(r"^\*\*`(?P<mod>[a-z_]+)\.(?P<name>[A-Za-z_]\w*)(?P<sig>\(.*\))?`\*\*"
                   r"(?P<rest>.*)$", re.M)


def _module(name):
    return importlib.import_module("isoform_dominance." + name)


def public(name):
    """``{name: object}`` of the module's public functions and classes."""
    m = _module(name)
    return {k: v for k, v in vars(m).items()
            if not k.startswith("_") and (inspect.isfunction(v) or inspect.isclass(v))
            and v.__module__ == m.__name__}


def _base(cls):
    if cls.__module__ == "builtins":
        return cls.__name__
    mod = cls.__module__.replace("isoform_dominance.", "")
    return "%s.%s" % (mod, cls.__name__)


def expected_line(mod, name, obj):
    if inspect.isclass(obj):
        return "**`%s.%s`** (a subclass of %s)" % (
            mod, name, " and ".join("`%s`" % _base(b) for b in obj.__bases__))
    return "**`%s.%s%s`**" % (mod, name, inspect.signature(obj))


def entries(text):
    """``[(module, name, the whole line)]`` of every entry line of the page."""
    return [(m["mod"], m["name"], m.group(0)) for m in ENTRY.finditer(text)]


def references(text):
    """``(module, name)`` of every ``module.name`` written in code on the page."""
    spans = re.findall(r"`([^`\n]+)`", text)
    pattern = re.compile(r"(?<![\w./:-])(%s)\.([A-Za-z_]\w*)" % "|".join(REFERABLE))
    return {(m.group(1), m.group(2)) for s in spans for m in pattern.finditer(s)}


TEXT = API.read_text(encoding="utf-8")
NAMES = [(mod, name) for mod in MODULES for name in sorted(public(mod))]


def test_the_count_is_the_reviewers_plus_this_round():
    # 98 functions and 11 classes, and io.groups_shape_error, which 2.9.0 adds
    kinds = [inspect.isclass(public(m)[n]) for m, n in NAMES]
    assert (kinds.count(False), kinds.count(True)) == (99, 11)


@pytest.mark.parametrize("mod,name", NAMES, ids=["%s.%s" % x for x in NAMES])
def test_each_public_name_has_one_entry_with_the_codes_signature(mod, name):
    lines = [line for m, n, line in entries(TEXT) if (m, n) == (mod, name)]
    assert lines == [expected_line(mod, name, public(mod)[name])]


def test_each_entry_is_under_its_modules_heading():
    where = None
    for line in TEXT.splitlines():
        head = re.match(r"^## `(\w+)`$", line)
        if head:
            where = head.group(1)
        m = ENTRY.match(line)
        if m:
            assert m["mod"] == where, line


def test_every_entry_and_reference_names_what_the_module_has():
    found = {(m, n) for m, n, _ in entries(TEXT)} | references(TEXT)
    assert len(found) > len(NAMES)
    absent = sorted("%s.%s" % (m, n) for m, n in found if not hasattr(_module(m), n))
    assert absent == []


def test_a_host_name_is_no_reference():
    text = ("`https://e<N>.rest.ensembl.org` `jul2023.rest.ensembl.org` "
            "`https://www.ensembl.org/help` `ensembl.resolve_server` `io.load_config(path)`")
    assert references(text) == {("ensembl", "resolve_server"), ("io", "load_config")}


def test_the_checks_can_fail():
    # an entry with an old signature, a missing entry and a name the module lacks
    stale = TEXT.replace(expected_line("ensembl", "resolve_server",
                                       public("ensembl")["resolve_server"]),
                         "**`ensembl.resolve_server(release=None, **retry)`**")
    assert [ln for m, n, ln in entries(stale) if (m, n) == ("ensembl", "resolve_server")] != [
        expected_line("ensembl", "resolve_server", public("ensembl")["resolve_server"])]
    assert ("annotate", "no_such_name") in references("`annotate.no_such_name`")
    assert not hasattr(_module("annotate"), "no_such_name")


def test_the_cli_is_named_as_internal():
    assert re.search(r"^## `cli`$", TEXT, re.M)
    section = TEXT[TEXT.index("## `cli`"):]
    assert "`isoform-dominance`" in section and "internal" in section
