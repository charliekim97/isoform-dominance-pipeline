"""A class made only of models whose 3' end was not found (GENCODE's ``mRNA_end_NF``).

``annotate`` clusters transcripts by their 3' terminal-exon acceptor.  The last annotated
exon of an ``mRNA_end_NF`` model is not known to be its 3' terminal exon -- it is often an
internal exon of a complete transcript -- so a class of such models alone can be a fragment
set against a complete class.  Only the GTF carries the tag: REST's lookup has no tags.  So
``annotate --gtf`` says so on stderr, leaving the config as it was, and
``identifiability --gtf`` gives every class ``incomplete_3prime_only``, null for sequence
that did not come from a GTF (REST's ``lookup/id/<id>?expand=1`` of an ``mRNA_end_NF``
transcript has no tag field).

On the GENCODE 50 extract (``tests/data/gencode_mini``): GSK3B's proposed comparison is
``iso_69aa`` against ``iso_400aa``, and the three transcripts of ``iso_69aa`` are all
``mRNA_end_NF`` while ``iso_400aa`` has none; one of STK11's twelve ``iso_433aa``
transcripts, ENST00000585851, is.
"""
import json
import os
import urllib.request

import pytest
from test_gtf_identifiability import FASTA, GTF, REST, Rest116
from test_gtf_parity import BY_NAME

from isoform_dominance import annotate, cli, io
from isoform_dominance import annotation_files as af
from isoform_dominance import identifiability as I

GOLDEN = os.path.join(os.path.dirname(__file__), "data", "v260_outputs")
GSK3B_69 = ["ENST00000678787", "ENST00000679066", "ENST00000679188"]
KEY = "incomplete_3prime_only"


@pytest.fixture
def rest(monkeypatch):
    """REST at release 116 from the recorded lookups, with their cDNA
    (``test_gtf_identifiability``'s), for ``identifiability``."""
    fake = Rest116()
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    return fake


@pytest.fixture
def annotate_rest(monkeypatch):
    """``annotate``'s REST calls answered from the recorded lookups, as
    ``test_pinned_outputs`` answers them."""
    def get(path, **kw):
        if path.startswith("/info/data"):
            return {"releases": [REST["release"]]}
        sym = path.split("?")[0].split("/")[4] if "/symbol/" in path else None
        if path.startswith("/lookup/symbol/"):
            return REST["lookups"][REST["lookup_symbol"].get(sym) or BY_NAME[sym][0]["id"]]
        if path.startswith("/xrefs/symbol/"):
            return [{"type": "gene", "id": g["id"]} for g in BY_NAME.get(sym, [])]
        if path.startswith("/lookup/id/"):
            return REST["lookups"][path.split("/")[3].split("?")[0]]
        raise AssertionError(path)
    monkeypatch.setattr(annotate, "_get", get)
    monkeypatch.setattr(annotate, "_post", lambda path, body, **kw:
                        {i: REST["lookups"][i] for i in body["ids"] if i in REST["lookups"]})
    monkeypatch.setattr(annotate.ensembl, "resolve_server", lambda release=None, **kw: None)


def _golden_bytes(name):
    with open(os.path.join(GOLDEN, name), "rb") as f:
        return f.read()


def _golden(name):
    """As a config file is written in text mode: in the platform's newline."""
    return _golden_bytes(name).replace(b"\n", os.linesep.encode())


def _annotate(tmp_path, capsys, gene, *extra):
    out = tmp_path / ("%s.json" % gene)
    rc = cli.main(["annotate", "--gene", gene, "--out", str(out)] + list(extra))
    got = capsys.readouterr()
    return rc, out, got.out, got.err


def _notes(err):
    return [line for line in err.splitlines() if "mRNA_end_NF" in line]


# --------------------------------------------------------------------------- #
# the record carries the tag; what annotate reads of it does not
# --------------------------------------------------------------------------- #
def test_the_gtf_record_carries_the_tags():
    gene = af.scan(GTF, symbol="GSK3B")[0]
    tagged = sorted(t["id"] for t in gene["Transcript"] if "mRNA_end_NF" in t["_tags"])
    assert set(GSK3B_69) <= set(tagged)
    canonical = next(t for t in gene["Transcript"] if t["is_canonical"])
    assert "Ensembl_canonical" in canonical["_tags"]
    # transcripts_of reads none of it: its records are REST's
    info = annotate.transcripts_of(gene, "GSK3B")
    assert {k for t in info["transcripts"] for k in t} == {
        "id", "protein_aa", "terminal_acceptor", "is_canonical"}


# --------------------------------------------------------------------------- #
# annotate --gtf: a NOTE on stderr, the config and --json as they were
# --------------------------------------------------------------------------- #
def test_annotate_gtf_notes_the_class_of_incomplete_models(tmp_path, capsys):
    rc, out, _, err = _annotate(tmp_path, capsys, "GSK3B", "--gtf", GTF)
    assert rc == 0
    assert out.read_bytes() == _golden("annotate_gtf_GSK3B.json")
    cfg = json.loads(out.read_text())
    assert cfg["primary_comparison"] == ["iso_69aa", "iso_400aa"]
    notes = _notes(err)
    assert len(notes) == 1, err
    note = notes[0]
    assert note.startswith("  NOTE: iso_69aa ")
    assert "iso_400aa" in note                      # the class it is compared with
    for t in GSK3B_69:
        assert t in note
    assert "3' terminal exon" in note


def test_annotate_gtf_json_is_the_config_as_it_was(tmp_path, capsys):
    rc, out, stdout, err = _annotate(tmp_path, capsys, "GSK3B", "--gtf", GTF, "--json")
    assert rc == 0
    assert stdout == _golden_bytes("annotate_gtf_GSK3B.json").decode() + "\n"
    assert len(_notes(err)) == 1


def test_annotate_from_rest_says_nothing_of_it(annotate_rest, tmp_path, capsys):
    """REST's lookup has no tags: nothing to say, and the config is REST's as it was."""
    rc, out, _, err = _annotate(tmp_path, capsys, "GSK3B")
    assert rc == 0
    assert out.read_bytes() == _golden("annotate_rest_GSK3B.json")
    assert _notes(err) == []
    rc, _, stdout, err = _annotate(tmp_path, capsys, "GSK3B", "--json")
    assert stdout == _golden_bytes("annotate_rest_GSK3B.json").decode() + "\n"
    assert _notes(err) == []


@pytest.mark.parametrize("gene", ["LEPR", "FOXO1", "STK11", "AXIN1", "CD99"])
def test_annotate_gtf_notes_no_class_with_a_complete_model(tmp_path, capsys, gene):
    """STK11's iso_433aa holds one mRNA_end_NF model of twelve: a complete one is there."""
    rc, out, _, err = _annotate(tmp_path, capsys, gene, "--gtf", GTF)
    assert rc == 0
    assert out.read_bytes() == _golden("annotate_gtf_%s.json" % gene)
    assert _notes(err) == []


# --------------------------------------------------------------------------- #
# identifiability --gtf: incomplete_3prime_only per class, null where no GTF said
# --------------------------------------------------------------------------- #
def _cfg(gene):
    return annotate.build_config_from_gtf(gene, GTF)


def test_identifiability_gtf_marks_the_class_of_incomplete_models():
    res = I.analyze(_cfg("GSK3B"), gtf=GTF, transcripts_fasta=FASTA)
    assert {g: r[KEY] for g, r in res["groups"].items()} == {"iso_69aa": True,
                                                            "iso_400aa": False}


def test_one_complete_model_makes_the_class_complete():
    res = I.analyze(_cfg("STK11"), gtf=GTF, transcripts_fasta=FASTA)
    assert {g: r[KEY] for g, r in res["groups"].items()} == {"iso_433aa_alt": False,
                                                            "iso_433aa": False}
    assert "ENST00000585851" in _cfg("STK11")["groups"]["iso_433aa"]


def test_every_class_has_the_key_and_the_note_is_of_the_compared_ones(tmp_path, capsys):
    """A third class, not compared, of incomplete models: its key is true and no NOTE is
    said of it; the compared classes are complete."""
    cfg = _cfg("GSK3B")
    cfg["groups"] = {"incomplete": cfg["groups"]["iso_69aa"][:1],
                     "a": cfg["groups"]["iso_400aa"][:6], "b": cfg["groups"]["iso_400aa"][6:]}
    cfg["primary_comparison"] = ["a", "b"]
    path = tmp_path / "cfg.json"
    io.dump_json(cfg, str(path))
    rc = cli.main(["identifiability", "--config", str(path), "--gtf", GTF,
                   "--transcripts-fasta", FASTA, "--json"])
    got = capsys.readouterr()
    assert rc in (0, 3)
    groups = json.loads(got.out)["groups"]
    assert {g: r[KEY] for g, r in groups.items()} == {"incomplete": True, "a": False,
                                                     "b": False}
    assert _notes(got.err) == []


def test_identifiability_gtf_notes_the_compared_class(tmp_path, capsys):
    path = tmp_path / "cfg.json"
    io.dump_json(_cfg("GSK3B"), str(path))
    rc = cli.main(["identifiability", "--config", str(path), "--gtf", GTF,
                   "--transcripts-fasta", FASTA])
    err = capsys.readouterr().err
    assert rc in (0, 3)
    notes = _notes(err)
    assert len(notes) == 1, err
    assert notes[0].startswith("  NOTE: iso_69aa: ")
    assert "iso_400aa" not in notes[0]
    for t in GSK3B_69:
        assert t in notes[0]


def test_the_verdict_and_figures_do_not_move():
    """The key is all the GTF adds: the report less it is the report from the same
    sequence supplied."""
    cfg = _cfg("GSK3B")
    by_files = I.analyze(cfg, gtf=GTF, transcripts_fasta=FASTA, min_log2fc=0.5)
    seqs = af.sequences_for(af.scan(GTF, symbol="GSK3B")[0], FASTA)
    needed = [t for ids in cfg["groups"].values() for t in ids]
    supplied = I.analyze(cfg, sequences={t: seqs[t] for t in needed},
                         background_sequences={t: s for t, s in seqs.items()
                                               if t not in needed},
                         background_gene_transcripts=False, min_log2fc=0.5)
    assert all(r[KEY] is None for r in supplied["groups"].values())
    for res in (by_files, supplied):
        for r in res["groups"].values():
            r.pop(KEY)
    assert json.dumps(by_files["groups"], sort_keys=True, default=str) \
        == json.dumps(supplied["groups"], sort_keys=True, default=str)
    assert (by_files["verdict"], by_files["contrast"]) == (supplied["verdict"],
                                                          supplied["contrast"])


def test_rest_and_supplied_sequence_say_null(rest):
    cfg = _cfg("GSK3B")
    res = I.analyze(cfg)
    assert rest.calls                               # REST, not the files
    assert {g: r[KEY] for g, r in res["groups"].items()} == {"iso_69aa": None,
                                                            "iso_400aa": None}


def test_the_fasta_alone_says_null():
    cfg = _cfg("GSK3B")
    res = I.analyze(cfg, transcripts_fasta=FASTA, background_gene_transcripts=False)
    assert {g: r[KEY] for g, r in res["groups"].items()} == {"iso_69aa": None,
                                                            "iso_400aa": None}


def test_a_rerun_from_saved_inputs_says_null(tmp_path, capsys):
    path, saved = tmp_path / "cfg.json", tmp_path / "inputs.json"
    io.dump_json(_cfg("GSK3B"), str(path))
    assert cli.main(["identifiability", "--config", str(path), "--gtf", GTF,
                     "--transcripts-fasta", FASTA, "--save-inputs", str(saved),
                     "--json"]) in (0, 3)
    first = json.loads(capsys.readouterr().out)["groups"]
    assert first["iso_69aa"][KEY] is True
    assert cli.main(["identifiability", "--config", str(path), "--inputs", str(saved),
                     "--json"]) in (0, 3)
    again = json.loads(capsys.readouterr().out)["groups"]
    assert {g: r[KEY] for g, r in again.items()} == {"iso_69aa": None, "iso_400aa": None}

