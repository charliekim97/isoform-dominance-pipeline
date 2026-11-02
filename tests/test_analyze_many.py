"""``identifiability.analyze_many``: many configs against one background FASTA.

``analyze_many(jobs, batch_size=16, errors="raise", **shared)`` gives, for each job, the
report ``analyze(config, **shared, **own)`` gives, in job order.  A batch of jobs reads the
FASTA once per pass for all of them -- a first pass, a second pass when a job needs one,
both again for each further round of the loop that takes a gene-background transcript's
sequence from the FASTA, and one header pass -- and takes its SHA-256 once.

Reports are compared as ``json.dumps``: a class with no informative read puts NaN in
``unique_read_counting_se``, and NaN is not equal to itself, so two equal reports can be
unequal by ``==``.
"""
import gzip
import json
import random
import urllib.request

import pytest

from isoform_dominance import annotate, index_scope, io
from isoform_dominance import identifiability as I
import test_decoys as DC
import test_fasta_competitors as FC
import test_identical_background as IB
import test_identical_copies as IC
import test_whole_index_rule as W
from test_gtf_identifiability import FASTA as MINI_FASTA
from test_gtf_identifiability import GTF, RECORDS as MINI_RECORDS


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("no request is made here")
    monkeypatch.setattr(urllib.request, "urlopen", refuse)


@pytest.fixture
def small_genome(monkeypatch):
    """A record longer than 20 kb counts as long here, as in test_decoys."""
    monkeypatch.setattr(I, "LONG_RECORD", 20_000)


class _Ids:
    """A fixture's transcript ids, each renamed to one of this module's own."""

    def __init__(self, n):
        self.n, self.map = n, {}

    def __call__(self, t):
        return self.map.setdefault(t, "ENST%02d%09d" % (self.n, len(self.map) + 1))

    def seqs(self, d):
        return {self(t): s for t, s in d.items()}

    def cfg(self, cfg, **kw):
        return dict(cfg, groups={g: [self(t) for t in ids] for g, ids in cfg["groups"].items()},
                    **kw)


def _gencode(tid, seq, gene="GENEX", gene_id="ENSG00000000077"):
    return ">%s.1|%s.1|-|-|%s-201|%s|%d|protein_coding|\n%s\n" % (
        tid, gene_id, gene, gene, len(seq), seq)


PLAIN = {"background_gene_transcripts": False}
CFG2 = {"gene": "GENEX", "groups": {"A": [IC.T1], "B": [IC.T2]},
        "primary_comparison": ["A", "B"]}


def _system(n, seed, records, gene="GENEX", gene_id="ENSG00000000077"):
    """test_whole_index_rule's gene, background and FASTA records, renamed: (config, own)."""
    ids = _Ids(n)
    cfg, seqs, bg, recs = W._system(random.Random(seed))
    for t, s in recs.items():
        records.append(_gencode(ids(t), s, gene, gene_id))
    return (ids.cfg(cfg, gene=gene),
            dict(PLAIN, sequences=ids.seqs(seqs), background_sequences=ids.seqs(bg))), ids


def _loop_job(records, n=20, seed=13):
    """A config whose gene-background transcript the FASTA holds with other sequence, and
    whose configured transcripts the FASTA holds under another gene name than the config's,
    beside a gene of that name under another id: a second round, the headers read twice by
    analyze, and a same-name copy."""
    (cfg, own), ids = _system(n, seed, records, gene="LOOPX", gene_id="ENSG00000000020")
    bg = own["background_sequences"]
    t = sorted(bg)[0]
    records.append(_gencode(t, bg[t][:50] + W._rand(random.Random(seed), 70), "LOOPX",
                            "ENSG00000000020"))
    for tid in sorted(own["sequences"]):
        records.append(_gencode(tid, own["sequences"][tid], "LOOPX", "ENSG00000000020"))
    records.append(_gencode("ENST%02d900000001" % n, W._rand(random.Random(seed + 1), 90),
                            "LOOPX", "ENSG00000000021"))
    return dict(cfg, gene="LOOPY"), own


def _mixed(tmp_path):
    """``(fasta, [(config, own)])``: the FASTA situations of the existing tests, in one FASTA."""
    records, jobs = [], []
    # test_whole_index_rule: paralogues, a repeat many records hold, copies, and records
    # that share an added record's own part (a second pass)
    for n, seed in enumerate((3, 5, 8), start=1):
        jobs.append(_system(n, seed, records)[0])
    jobs.append((jobs[0][0], dict(jobs[0][1], estimand="abundance", min_log2fc=0.5)))
    jobs.append((jobs[1][0], dict(jobs[1][1], keep_duplicates=True)))
    jobs.append((jobs[2][0], dict(jobs[2][1], max_window_records=1, min_log2fc=0.5,
                                  estimand="abundance")))
    jobs.append(_loop_job(records))
    # test_identical_copies: a record with a configured transcript's sequence
    ids = _Ids(5)
    records.append(_gencode(ids(IC.TWIN), IC.SEQS[IC.T1]))
    jobs.append((ids.cfg(CFG2), dict(PLAIN, sequences=ids.seqs(IC.SEQS))))
    # test_identical_background: a gene-background twin, and a record with the sequence of
    # a background transcript
    ids = _Ids(6)
    records.append(_gencode(ids("ENST00000000099"), IB.BG[IB.OTHER]))
    jobs.append((ids.cfg(IB.CFG), dict(PLAIN, sequences=ids.seqs(IB.SEQS),
                                       background_sequences=ids.seqs(IB.BG), min_log2fc=0.5)))
    # test_fasta_competitors: a competitor, and a record that shares nothing
    ids = _Ids(7)
    records.append(_gencode(ids(FC.COMP), FC.COMPETITOR))
    records.append(_gencode(ids(FC.STRANGER), FC._seq(900, 5)))
    jobs.append((ids.cfg(FC.CFG), dict(PLAIN, sequences=ids.seqs(FC.SEQS), min_log2fc=0.5)))
    # test_decoys: the gene's transcripts and a chromosome that holds its exons, no decoys
    ids = _Ids(8)
    records.append(_gencode(ids("ENST00000000009"), DC.OTHER))
    records.append(">chrF\n%s\n" % DC.CHROM)
    jobs.append((ids.cfg(DC.CFG), dict(PLAIN, sequences=ids.seqs(DC.SEQS))))
    # every column shorter than the window: no first pass for it
    ids = _Ids(9)
    jobs.append((ids.cfg(CFG2), dict(PLAIN, sequences=ids.seqs(
        {IC.T1: IC.SEQS[IC.T1][:20], IC.T2: IC.SEQS[IC.T2][:25]}))))
    # local annotation: GENCODE 50's LEPR and CD99, whose records are in the FASTA too
    records.extend(">%s\n%s\n" % (h, s) for h, s in MINI_RECORDS)
    for gene in ("LEPR", "CD99"):
        jobs.append((annotate.build_config_from_gtf(gene, GTF),
                     {"gtf": GTF, "transcripts_fasta": MINI_FASTA, "min_log2fc": 0.5}))
    path = tmp_path / "mixed.fa"
    path.write_text("".join(records))
    return str(path), jobs


def _with_inputs_out(jobs):
    """The jobs, each with an ``inputs_out`` dict of its own, and those dicts."""
    outs = [{} for _ in jobs]
    return [(c, dict(o, inputs_out=d)) for (c, o), d in zip(jobs, outs, strict=True)], outs


def _same_as_analyze(jobs, got, outs, shared):
    for (cfg, own), mine, mine_out in zip(jobs, got, outs, strict=True):
        want_out = {}
        want = I.analyze(cfg, **shared, **dict(own, inputs_out=want_out))
        assert json.dumps(mine) == json.dumps(want)
        assert json.dumps(mine_out) == json.dumps(want_out)


# --------------------------------------------------------------------------- #
# each job, analyze's report
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("batch_size", [32, 5])
def test_each_job_gets_the_report_analyze_gives(tmp_path, small_genome, batch_size):
    fa, jobs = _mixed(tmp_path)
    shared = {"background_fasta": fa}
    with_out, outs = _with_inputs_out(jobs)
    got = list(I.analyze_many(with_out, batch_size=batch_size, **shared))
    assert len(got) == len(jobs)
    _same_as_analyze(with_out, got, outs, shared)
    # and the situations are there: nothing above is compared on empty fields
    bgs = [r["background"] for r in got]
    assert any(b["sequence_from_fasta"] for b in bgs)
    assert any(b["fasta_competitors"] for b in bgs)
    assert any(b["fasta_left_out"] for b in bgs)
    assert any(b["windows_dropped"]["total"] for b in bgs)
    assert any(b["fasta_long_records"] for b in bgs)
    assert any(b["same_name_copies"] for b in bgs)
    assert any(src == "fasta" for b in bgs for src in (b["identical_source"] or {}).values())
    assert any(r["annotation"]["source"] for r in got)
    assert {r["decision_estimand"] for r in got} == {"fragment", "abundance"}
    assert {r["effect_resolvable"] is None for r in got} == {True, False}
    assert any(not b["windows_dropped"]["configured_of"] for b in bgs)


def test_with_decoys_each_job_gets_the_report_analyze_gives(tmp_path, small_genome):
    fa, jobs = tmp_path / "gentrome.fa", []
    records = [DC.GENTROME]
    for n in (1, 2):
        jobs.append(_system(n, 20 + n, records)[0])
    jobs.append((DC.CFG, dict(PLAIN, sequences=DC.SEQS, min_log2fc=0.5)))
    fa.write_text("".join(records))
    shared = {"background_fasta": str(fa), "decoys": DC._decoys(tmp_path, "chrF", "chrG",
                                                                 "chrAbsent")}
    with_out, outs = _with_inputs_out(jobs)
    got = list(I.analyze_many(with_out, **shared))
    _same_as_analyze(with_out, got, outs, shared)
    assert {r["background"]["decoys_skipped"] for r in got} == {2}
    assert {tuple(r["background"]["decoys_absent"]) for r in got} == {("chrAbsent",)}


def test_without_a_fasta_it_is_analyze_one_by_one():
    jobs = [(CFG2, dict(PLAIN, sequences=IC.SEQS)),
            (IB.CFG, dict(PLAIN, sequences=IB.SEQS, background_sequences=IB.BG,
                          estimand="abundance", min_log2fc=0.5))]
    with_out, outs = _with_inputs_out(jobs)
    got = list(I.analyze_many(with_out, batch_size=1))
    _same_as_analyze(with_out, got, outs, {})


def test_a_config_alone_is_a_job_and_shared_arguments_reach_every_job(tmp_path):
    fa = tmp_path / "bg.fa"
    fa.write_text(_gencode(IC.TWIN, IC.SEQS[IC.T1]))
    shared = {"sequences": IC.SEQS, "background_gene_transcripts": False,
              "background_fasta": str(fa), "tpm": 3.0}
    got = list(I.analyze_many([CFG2, (CFG2, {"min_log2fc": 0.4}), [CFG2, {}]], **shared))
    assert [json.dumps(r) for r in got] == [
        json.dumps(I.analyze(CFG2, **shared)), json.dumps(I.analyze(CFG2, min_log2fc=0.4,
                                                                     **shared)),
        json.dumps(I.analyze(CFG2, **shared))]


# --------------------------------------------------------------------------- #
# how often a batch reads the FASTA
# --------------------------------------------------------------------------- #
class _Reads:
    """What reads the FASTA at ``path``: the opens through ``io.open_text`` (the whole-index
    passes) and through ``index_scope.open_text`` (the header pass), and the calls of
    ``io.file_sha256``.  Other files are not counted."""

    def __init__(self, monkeypatch, path):
        self.path, self.passes, self.headers, self.sha = str(path), 0, 0, 0
        real_open, real_sha = io.open_text, io.file_sha256

        def opener(what):
            def call(p, *a, **k):
                if str(p) == self.path:
                    setattr(self, what, getattr(self, what) + 1)
                return real_open(p, *a, **k)
            return call

        def sha(p, *a, **k):
            if str(p) == self.path:
                self.sha += 1
            return real_sha(p, *a, **k)

        monkeypatch.setattr(io, "open_text", opener("passes"))
        monkeypatch.setattr(index_scope, "open_text", opener("headers"))
        monkeypatch.setattr(io, "file_sha256", sha)

    def now(self):
        return self.passes, self.headers, self.sha


def _alone(monkeypatch, reads, job, shared):
    """analyze's reads for one job: ``[(pass 1, pass 2)]`` per round, header passes, SHA-256."""
    cfg, own = job
    rounds, real = [], I._whole_index

    def spy(*args, **kwargs):
        before = reads.passes
        out = real(*args, **kwargs)
        first = int(bool(set(args[1])))
        rounds.append((first, reads.passes - before - first))
        return out

    monkeypatch.setattr(I, "_whole_index", spy)
    _, h, s = reads.now()
    I.analyze(cfg, **shared, **own)
    monkeypatch.setattr(I, "_whole_index", real)
    return rounds, reads.headers - h, reads.sha - s


def _contract(alone):
    """The reads of a batch of these jobs, by the contract: per round a first pass when a
    job of the round has a column window, and a second when one needs it; one header pass
    and one SHA-256."""
    n = max(len(r) for r, _, _ in alone)
    passes = 0
    for i in range(n):
        passes += max(r[i][0] for r, _, _ in alone if len(r) > i)
        passes += max(r[i][1] for r, _, _ in alone if len(r) > i)
    return passes, 1, 1


def _batch_reads(reads, jobs, shared, **kw):
    before = reads.now()
    got = list(I.analyze_many(jobs, **kw, **shared))
    return tuple(a - b for a, b in zip(reads.now(), before, strict=True)), got


def _read_cases(tmp_path):
    records = []
    loop = _loop_job(records)
    second = _system(1, 3, records)[0]
    others = [_system(n, seed, records)[0] for n, seed in ((2, 5), (3, 8), (4, 11), (5, 2))]
    ids = _Ids(9)
    records.append(_gencode(ids(IC.TWIN), IC.SEQS[IC.T1]))
    others += [(ids.cfg(CFG2), dict(PLAIN, sequences=ids.seqs(IC.SEQS)))]
    ids = _Ids(10)
    others += [(ids.cfg(CFG2), dict(PLAIN, sequences=ids.seqs(
        {IC.T1: IC.SEQS[IC.T1][:20], IC.T2: IC.SEQS[IC.T2][:25]})))]
    fa = tmp_path / "reads.fa"
    fa.write_text("".join(records))
    return str(fa), loop, second, others


def test_a_batch_reads_the_fasta_as_often_with_two_jobs_as_with_twelve(tmp_path, monkeypatch):
    fa, loop, second, others = _read_cases(tmp_path)
    shared = {"background_fasta": fa}
    reads = _Reads(monkeypatch, fa)
    small = [loop, second]
    big = small + others + [(c, dict(o, min_log2fc=0.5)) for c, o in others[:4]]
    assert len(big) == 12 and big[7] == others[-1]
    alone = [_alone(monkeypatch, reads, j, shared) for j in big]
    # analyze reads per config as it did: a first pass when there is a column window, a
    # second when needed, both again for each further round, the headers once or twice
    loop_alone, second_alone, short_alone = alone[0], alone[1], alone[7]
    assert [first for first, _ in loop_alone[0]] == [1, 1] and loop_alone[2] == 1
    assert loop_alone[1] == 2              # the header gives its transcripts another name
    assert second_alone[0][0] == (1, 1)
    assert short_alone == ([(0, 0)], 1, 1)
    assert all(len(r) == 1 and h in (1, 2) and s == 1 for r, h, s in alone[2:])
    n_alone = sum(sum(a + b for a, b in r) + h for r, h, _ in alone)

    want = _contract(alone[:2])
    assert _contract(alone) == want
    got_small, _ = _batch_reads(reads, small, shared)
    got_big, reports = _batch_reads(reads, big, shared)
    assert got_small == got_big == want
    assert sum(want[:2]) < n_alone
    assert [json.dumps(r) for r in reports] == [
        json.dumps(I.analyze(c, **shared, **o)) for c, o in big]


def test_a_batch_of_jobs_that_need_no_second_round_reads_each_pass_once(tmp_path,
                                                                        monkeypatch):
    fa, _, second, others = _read_cases(tmp_path)
    shared = {"background_fasta": fa}
    reads = _Reads(monkeypatch, fa)
    jobs = [second] + others
    alone = [_alone(monkeypatch, reads, j, shared) for j in jobs]
    assert _contract(alone) == (2, 1, 1)
    assert _batch_reads(reads, jobs, shared)[0] == (2, 1, 1)


def test_a_batch_whose_jobs_have_no_column_window_reads_only_the_headers(tmp_path,
                                                                        monkeypatch):
    fa, _, _, others = _read_cases(tmp_path)
    short = others[-1]
    reads = _Reads(monkeypatch, fa)
    shared = {"background_fasta": fa}
    assert _batch_reads(reads, [short, short, short], shared)[0] == (0, 1, 1)


def test_the_batches_are_as_many_as_batch_size_makes_them(tmp_path, monkeypatch):
    fa, _, second, others = _read_cases(tmp_path)
    shared = {"background_fasta": fa}
    jobs = ([second] + others)[:7]
    reads = _Reads(monkeypatch, fa)
    for size, batches in ((3, 3), (7, 1), (1, 7), (16, 1)):
        got, _ = _batch_reads(reads, jobs, shared, batch_size=size)
        assert got[1:] == (batches, batches)


def test_a_batch_computes_one_job_at_a_time(tmp_path, monkeypatch):
    fa, _, second, others = _read_cases(tmp_path)
    shared = {"background_fasta": fa}
    order = []
    real_matrix, real_sha = I.compatibility_matrix, io.file_sha256

    def matrix(*a, **k):
        order.append("matrix")
        return real_matrix(*a, **k)

    def sha(p, *a, **k):
        order.append("sha256")
        return real_sha(p, *a, **k)

    monkeypatch.setattr(I, "compatibility_matrix", matrix)
    monkeypatch.setattr(io, "file_sha256", sha)
    list(I.analyze_many([second, others[0], others[1]], **shared))
    # the SHA-256 is answered where the first job asks, not held for the batch: the first
    # job ends before the second builds its system
    assert order == ["matrix", "sha256", "matrix", "matrix"]


# --------------------------------------------------------------------------- #
# the contract: shapes, arguments, errors
# --------------------------------------------------------------------------- #
def _small(tmp_path):
    fa = tmp_path / "bg.fa"
    fa.write_text(_gencode(IC.TWIN, IC.SEQS[IC.T1]) + _gencode(FC.COMP, FC.COMPETITOR))
    return str(fa), {"sequences": dict(IC.SEQS), "background_gene_transcripts": False}


def _outcome(call, *a, **k):
    try:
        return "report", json.dumps(call(*a, **k))
    except Exception as e:
        return type(e), e.args


def _yielded(result):
    if isinstance(result, BaseException):
        return type(result), result.args
    return "report", json.dumps(result)


BAD_CONFIG = {"gene": "GENEX", "groups": {"A": [IC.T1]}, "primary_comparison": ["A", "B"]}


@pytest.mark.parametrize("job", [
    "ENST1", 5, None, (CFG2,), (CFG2, {}, {}), (CFG2, "own"), [CFG2, None], {CFG2["gene"]},
], ids=["str", "int", "none", "one", "three", "own-str", "own-none", "set"])
def test_a_job_of_another_shape_is_that_jobs_type_error(tmp_path, job):
    fa, kw = _small(tmp_path)
    got = list(I.analyze_many([CFG2, job, CFG2], errors="return", background_fasta=fa, **kw))
    assert type(got[1]) is TypeError and "config" in str(got[1])
    assert _yielded(got[0]) == _yielded(got[2]) == _outcome(I.analyze, CFG2,
                                                           background_fasta=fa, **kw)


@pytest.mark.parametrize("name", ["background_fasta", "decoys", "k", "window", "canonical"])
def test_an_argument_that_reads_the_fasta_is_shared_only(tmp_path, name):
    fa, kw = _small(tmp_path)
    got = list(I.analyze_many([CFG2, (CFG2, {name: None})], errors="return",
                              background_fasta=fa, **kw))
    assert type(got[1]) is ValueError and name in str(got[1])
    assert _yielded(got[0])[0] == "report"


def test_a_name_in_shared_and_own_is_the_type_error_analyze_raises(tmp_path):
    fa, kw = _small(tmp_path)
    got = list(I.analyze_many([(CFG2, {"tpm": 2.0}), (CFG2, {"no_such": 1})],
                              errors="return", background_fasta=fa, tpm=1.0, **kw))
    for mine, own in zip(got, ({"tpm": 2.0}, {"no_such": 1}), strict=True):
        with pytest.raises(TypeError):
            I.analyze(CFG2, background_fasta=fa, tpm=1.0, **kw, **own)
        assert type(mine) is TypeError


def test_the_checks_go_shape_then_fasta_arguments_then_names(tmp_path):
    fa, kw = _small(tmp_path)
    got = list(I.analyze_many([(CFG2, {"window": 31, "tpm": 1.0}), (CFG2, "x")],
                              errors="return", background_fasta=fa, tpm=1.0, **kw))
    assert [type(e) for e in got] == [ValueError, TypeError]


def test_inputs_out_is_own_only_and_said_at_the_first_next(tmp_path):
    fa, kw = _small(tmp_path)
    many = I.analyze_many([CFG2], background_fasta=fa, inputs_out={}, **kw)
    with pytest.raises(ValueError, match="inputs_out"):
        next(many)


@pytest.mark.parametrize("bad", [{"batch_size": 0}, {"batch_size": -1}, {"batch_size": True},
                                 {"batch_size": 2.0}, {"batch_size": "2"},
                                 {"errors": "ignore"}, {"errors": None}])
def test_a_bad_batch_size_or_errors_is_said_at_the_first_next(bad):
    many = I.analyze_many([CFG2], **bad)
    with pytest.raises(ValueError):
        next(many)


def _jobs_with_errors():
    return [(CFG2, {}), (BAD_CONFIG, {}), (CFG2, {"min_log2fc": 0.5}),
            (dict(CFG2, gene=["GENEX"]), {}), (CFG2, {"tpm": 2.0})]


@pytest.mark.parametrize("batch_size", [1, 2, 16])
def test_errors_return_yields_each_jobs_error_in_its_place(tmp_path, batch_size):
    fa, kw = _small(tmp_path)
    got = list(I.analyze_many(_jobs_with_errors(), batch_size=batch_size, errors="return",
                              background_fasta=fa, **kw))
    want = [_outcome(I.analyze, c, background_fasta=fa, **kw, **o)
            for c, o in _jobs_with_errors()]
    assert [_yielded(r) for r in got] == want
    # a config whose gene is a list fails in the header pass: in the joint call, for it alone
    assert want[3][0] is TypeError and want[1][0] is ValueError


@pytest.mark.parametrize("batch_size", [1, 2, 16])
def test_errors_raise_raises_after_the_reports_before_it(tmp_path, batch_size):
    fa, kw = _small(tmp_path)
    jobs = _jobs_with_errors()
    many = I.analyze_many(jobs, batch_size=batch_size, background_fasta=fa, **kw)
    assert _yielded(next(many)) == _outcome(I.analyze, CFG2, background_fasta=fa, **kw)
    with pytest.raises(ValueError) as e:
        next(many)
    assert e.value.args == _outcome(I.analyze, BAD_CONFIG, background_fasta=fa, **kw)[1]
    with pytest.raises(StopIteration):
        next(many)


def test_errors_raise_on_a_job_of_the_header_pass(tmp_path):
    fa, kw = _small(tmp_path)
    jobs = [(CFG2, {}), (dict(CFG2, gene=["GENEX"]), {}), (CFG2, {})]
    many = I.analyze_many(jobs, background_fasta=fa, **kw)
    assert _yielded(next(many)) == _outcome(I.analyze, CFG2, background_fasta=fa, **kw)
    with pytest.raises(TypeError) as e:
        next(many)
    with pytest.raises(TypeError) as want:
        I.analyze(dict(CFG2, gene=["GENEX"]), background_fasta=fa, **kw)
    assert e.value.args == want.value.args


def test_errors_raise_may_have_computed_the_jobs_after_it_in_its_batch(tmp_path):
    fa, kw = _small(tmp_path)
    after = {}
    many = I.analyze_many([(BAD_CONFIG, {}), (CFG2, {"inputs_out": after})], batch_size=2,
                          background_fasta=fa, **kw)
    with pytest.raises(ValueError):
        next(many)
    assert after["sequences"] == IC.SEQS                # computed, and not yielded


def test_an_interrupt_is_raised_whatever_errors_says(tmp_path, monkeypatch):
    fa, kw = _small(tmp_path)
    real, n = I.compatibility_matrix, []

    def matrix(*a, **k):
        n.append(1)
        if len(n) == 2:
            raise KeyboardInterrupt
        return real(*a, **k)

    monkeypatch.setattr(I, "compatibility_matrix", matrix)
    with pytest.raises(KeyboardInterrupt):
        list(I.analyze_many([CFG2, CFG2, CFG2], errors="return", background_fasta=fa, **kw))


def test_a_fasta_that_cannot_be_read_is_the_batchs_error_whatever_errors_says(tmp_path):
    _, kw = _small(tmp_path)
    many = I.analyze_many([CFG2, BAD_CONFIG], errors="return",
                          background_fasta=str(tmp_path / "absent.fa"), **kw)
    with pytest.raises(FileNotFoundError):
        next(many)


def test_jobs_are_pulled_a_batch_at_a_time(tmp_path):
    fa, kw = _small(tmp_path)
    pulled = []

    def jobs():
        for i in range(5):
            pulled.append(i)
            yield (CFG2, {"tpm": 1.0 + i})

    many = I.analyze_many(jobs(), batch_size=2, background_fasta=fa, **kw)
    assert pulled == []
    got = [next(many)]
    assert pulled == [0, 1]
    got.append(next(many))
    assert pulled == [0, 1]
    got += list(many)
    assert pulled == [0, 1, 2, 3, 4]
    assert [_yielded(r) for r in got] == [
        _outcome(I.analyze, CFG2, background_fasta=fa, tpm=1.0 + i, **kw) for i in range(5)]


def test_no_jobs_no_reports():
    assert list(I.analyze_many([])) == []
    assert list(I.analyze_many(iter(()), batch_size=3)) == []


def test_an_error_pulling_the_jobs_comes_after_the_reports_of_those_pulled(tmp_path):
    fa, kw = _small(tmp_path)

    def jobs():
        yield CFG2
        yield (CFG2, {"tpm": 2.0})
        yield (CFG2, {"tpm": 3.0})
        raise OSError("the job list broke")

    many = I.analyze_many(jobs(), batch_size=2, errors="return", background_fasta=fa, **kw)
    got = [next(many), next(many), next(many)]
    with pytest.raises(OSError, match="the job list broke"):
        next(many)
    assert [_yielded(r) for r in got] == [
        _outcome(I.analyze, CFG2, background_fasta=fa, tpm=t, **kw) for t in (10.0, 2.0, 3.0)]


# --------------------------------------------------------------------------- #
# what an error kept by the caller holds
# --------------------------------------------------------------------------- #
#: the tables of the joint calls: the header pass's, and each whole-index pass's
_TABLES = {"records", "by_transcript", "by_name", "by_first", "live", "plan", "owner", "find",
           "by_tid", "by_len"}


def _frames(error):
    tb = error.__traceback__
    while tb is not None:
        yield tb.tb_frame
        tb = tb.tb_next


def test_a_kept_error_holds_none_of_its_joint_calls_tables(tmp_path, monkeypatch):
    # a job's error in a joint call is yielded with its traceback, which holds the joint
    # call's frame for as long as the caller keeps the error: the frame must not hold the
    # batch's tables -- every header of the FASTA, or every job's state
    fa, kw = _small(tmp_path)
    real = I._WholeIndexJob.pass1

    def pass1(self, *a, **k):
        if self.max_window_records == 7:
            raise RuntimeError("this job's own")
        return real(self, *a, **k)

    monkeypatch.setattr(I._WholeIndexJob, "pass1", pass1)
    got = list(I.analyze_many([CFG2, (dict(CFG2, gene=["GENEX"]), {}),
                               (CFG2, {"max_window_records": 7}), CFG2],
                              errors="return", background_fasta=fa, **kw))
    assert [type(r) for r in got] == [dict, TypeError, RuntimeError, dict]
    joint = [[f for f in _frames(e) if f.f_code.co_name == name]
             for e, name in ((got[1], "_fasta_copies_many"), (got[2], "_whole_index_many"))]
    assert [len(frames) for frames in joint] == [1, 1]
    held = [sorted(n for n, v in f.f_locals.items() if n in _TABLES and v is not None)
            for frames in joint for f in frames]
    assert held == [[], []]


# --------------------------------------------------------------------------- #
# another window, strand or compression
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("other", [{"k": 25}, {"window": 25}, {"canonical": False},
                                   {"gzip": True}],
                         ids=["k25", "window25", "strand-aware", "gzipped"])
def test_each_job_gets_the_report_analyze_gives_at_any_window_strand_or_compression(
        tmp_path, other):
    other = dict(other)
    records = []
    jobs = [_system(n, seed, records)[0] for n, seed in ((1, 3), (2, 5), (3, 8))]
    jobs.append(_loop_job(records))
    text = "".join(records)
    if other.pop("gzip", False):
        path = tmp_path / "bg.fa.gz"
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write(text)
    else:
        path = tmp_path / "bg.fa"
        path.write_text(text)
    shared = dict(other, background_fasta=str(path))
    with_out, outs = _with_inputs_out(jobs)
    got = list(I.analyze_many(with_out, **shared))
    _same_as_analyze(with_out, got, outs, shared)
    assert any(r["background"]["fasta_competitors"] for r in got)


# --------------------------------------------------------------------------- #
# requests, and a GTF that changes while it is read
# --------------------------------------------------------------------------- #
def test_a_request_of_another_kind_is_refused_not_answered(monkeypatch):
    with pytest.raises(ValueError, match="kind"):
        I._answer(("another", ("x",), {}))

    def analysis(config, **kw):
        yield ("another", ("x",), {})

    monkeypatch.setattr(I, "_analysis", analysis)
    with pytest.raises(ValueError, match="kind"):
        list(I.analyze_many([CFG2], errors="return"))


def test_a_gtf_that_loses_the_gene_between_two_reads_is_a_value_error(monkeypatch):
    # the gene is found by its id under another name, so the GTF is read again by name; a
    # file that no longer gives the gene then is named, where 2.10.0 let StopIteration out
    cfg = dict(annotate.build_config_from_gtf("LEPR", GTF), gene="NOT_LEPR")
    files = {"gtf": GTF, "transcripts_fasta": MINI_FASTA}
    real, calls = I.af.scan, []

    def scan(*a, **k):
        calls.append(k)
        return real(*a, **k) if len(calls) % 2 == 1 else []

    monkeypatch.setattr(I.af, "scan", scan)
    with pytest.raises(ValueError, match="changed") as alone:
        I.analyze(cfg, **files)
    got = list(I.analyze_many([(cfg, files)], errors="return"))
    assert [(type(e), e.args) for e in got] == [(ValueError, alone.value.args)]
    assert len(calls) == 4                  # by id, then by name; for each call
