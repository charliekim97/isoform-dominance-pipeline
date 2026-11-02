"""``analyze`` runs ``_analysis``, a generator that asks for each read of a background FASTA.

The body of ``analyze`` is the generator ``_analysis``: where it needs the FASTA it yields a
request -- the whole-index passes, the header pass for same-name copies, the SHA-256 of the
file -- and is sent the answer.  ``analyze`` answers each request by calling the function it
always called, looked up on its module when the request comes, so that one config reads the
FASTA as it did through 2.10.0.  ``analyze_many`` answers the requests of several configs
together.
"""
import inspect
import json
import random

from isoform_dominance import identifiability as I
from isoform_dominance import index_scope, io


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


SHARED, ALT_A, ALT_B, GENE_OWN, OTHER = (_seq(600, 1), _seq(300, 2), _seq(300, 3),
                                         _seq(250, 4), _seq(400, 5))
T1, T2, BG = "ENST00000000001", "ENST00000000002", "ENST00000000011"
SEQS = {T1: SHARED + ALT_A, T2: SHARED + ALT_B}
CFG = {"gene": "GENEX", "groups": {"A": [T1], "B": [T2]}, "primary_comparison": ["A", "B"]}


def _fasta(tmp_path, records):
    p = tmp_path / "bg.fa"
    p.write_text("".join(">%s.1|ENSG00000000010.1|-|-|GENEX-201|GENEX|%d|protein_coding|\n%s\n"
                         % (t, len(s), s) for t, s in records.items()))
    return str(p)


def test_the_analysis_takes_the_arguments_analyze_takes():
    assert inspect.isgeneratorfunction(I._analysis)
    assert inspect.signature(I._analysis) == inspect.signature(I.analyze)


def test_analyze_hands_every_argument_to_the_analysis(monkeypatch):
    sig = inspect.signature(I.analyze)
    sentinels = {name: object() for name in sig.parameters}
    got = {}

    def fake(*args, **kwargs):
        got.update(sig.bind(*args, **kwargs).arguments)
        return "report"
        yield                                       # a generator, as _analysis is

    monkeypatch.setattr(I, "_analysis", fake)
    pos = [sentinels[n] for n, p in sig.parameters.items()
           if p.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD]
    kw = {n: sentinels[n] for n, p in sig.parameters.items()
          if p.kind is inspect.Parameter.KEYWORD_ONLY}
    assert I.analyze(*pos, **kw) == "report"
    assert set(got) == set(sentinels)
    assert all(got[n] is sentinels[n] for n in sentinels)


def test_the_analysis_asks_for_each_read_of_the_fasta_in_order(tmp_path):
    # the FASTA holds the gene-background transcript with other sequence: a second round
    fa = _fasta(tmp_path, {BG: GENE_OWN + SHARED[:200], "ENST00000000099": OTHER})
    kw = {"sequences": SEQS, "background_sequences": {BG: GENE_OWN},
          "background_gene_transcripts": False, "background_fasta": fa}
    gen = I._analysis(CFG, **kw)
    kinds, answer = [], None
    while True:
        try:
            kind, args, kwargs = gen.send(answer)
        except StopIteration as done:
            report = done.value
            break
        kinds.append(kind)
        assert args[0] == fa
        answer = {"whole_index": I._whole_index, "copies": index_scope.fasta_copies,
                  "sha256": io.file_sha256}[kind](*args, **kwargs)
    assert kinds == ["whole_index", "whole_index", "copies", "sha256"]
    assert json.dumps(report) == json.dumps(I.analyze(CFG, **kw))
    assert report["background"]["sequence_from_fasta"] == [BG]


def test_analyze_answers_with_the_function_its_module_has_when_asked(tmp_path, monkeypatch):
    fa = _fasta(tmp_path, {"ENST00000000099": SHARED[:300] + OTHER})
    calls = []

    def spy(name, real):
        def call(*args, **kwargs):
            calls.append((name, args[0], sorted(kwargs)))
            return real(*args, **kwargs)
        return call

    monkeypatch.setattr(I, "_whole_index", spy("whole", I._whole_index))
    monkeypatch.setattr(index_scope, "fasta_copies", spy("copies", index_scope.fasta_copies))
    monkeypatch.setattr(io, "file_sha256", spy("sha256", io.file_sha256))
    I.analyze(CFG, sequences=SEQS, background_gene_transcripts=False, background_fasta=fa)
    assert calls == [
        ("whole", fa, ["canonical", "decoys", "exclude_ids", "keep_duplicates",
                       "max_window_records"]),
        ("copies", fa, ["versions"]),
        ("sha256", fa, [])]
