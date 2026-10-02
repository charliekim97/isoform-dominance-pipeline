"""The ``--background-fasta`` scan gives 2.4.1's answer, faster.

2.4.1 folded every window of every record to its canonical form and looked it up, at about
2 Mb of sequence a second: GENCODE 50's reference-chromosome FASTA (644,292 records,
1,443,754,799 nt) took about nine minutes a run.  The scan now looks each query window up in
both orientations and reads a record's window only where a seed of it matches.  It is meant
to give the same answer, so it is tested against 2.4.1's own function, kept below as the
reference, on random cases that mix window lengths, both k-mer conventions, N, lower case,
windows planted in either orientation, records shorter than the window, identical copies,
excluded ids, and gzip or plain files.
"""
import gzip
import random

import pytest

from isoform_dominance import identifiability as I
from isoform_dominance import io


def _scan_241(path, query_kmers, k, canonical=True, exclude_ids=(),
              identical=None, identical_out=None):
    """``scan_background_fasta`` as released in 2.4.1, unchanged."""
    query = set(query_kmers)
    if not query:
        return set()
    exclude = {str(i).split(".")[0] for i in exclude_ids}
    identical = identical or {}
    lengths = {len(s) for s in identical}
    seen = set()

    def _consume(chunks, tid):
        if not chunks or tid in exclude:
            return
        seq = "".join(chunks).upper()
        if len(seq) in lengths and seq in identical:
            if identical_out is not None:
                identical_out[tid] = identical[seq]
            return
        for i in range(len(seq) - k + 1):
            km = seq[i:i + k]
            if canonical:
                km = I.canonical_kmer(km)
            if km in query:
                seen.add(km)

    chunks, tid = [], None
    with io.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                _consume(chunks, tid)
                head = line[1:].strip()
                tid = head.split("|")[0].split()[0].split(".")[0] if head else None
                chunks = []
            else:
                chunks.append(line.strip())
        _consume(chunks, tid)
    return seen


def _rand(r, n, alphabet="ACGT"):
    return "".join(r.choice(alphabet) for _ in range(n))


def _case(r, tmp_path, i):
    """One random case: the arguments both scans are called with."""
    k = r.choice([5, 6, 7, 11, 16, 17, 31, 32, 50, 101, r.randint(5, 101)])
    canonical = r.random() < 0.7
    seed = r.choice([1, 2, 3, 7, 16, 16, 16, k, k + 3, r.randint(1, k)])
    configured = {}
    for j in range(r.randint(1, 3)):
        body = _rand(r, r.randint(k, k + 250))
        if r.random() < 0.3:                     # an N, or a low-complexity stretch
            at = r.randrange(len(body))
            body = body[:at] + r.choice(["N", "A" * 40, "CA" * 30, "NNNN"]) + body[at:]
        configured["ENST%011d" % (j + 1)] = body
    query = set()
    for s in configured.values():
        query |= I.kmers(s, k, canonical)
    if r.random() < 0.2:                         # entries no record window can match
        extra = _rand(r, k)
        query.add(I.revcomp(extra) if canonical and I.canonical_kmer(extra) == extra
                  else extra)
        query.add(_rand(r, k + 1))

    def planted():
        src = r.choice(list(configured.values()))
        if len(src) < k:
            return ""
        a = r.randint(0, len(src) - k)
        w = src[a:a + r.randint(k, min(len(src) - a, 3 * k))]
        if r.random() < 0.5:
            w = I.revcomp(w)
        if r.random() < 0.2:
            w = w.lower()
        if r.random() < 0.15 and len(w) > 2:     # a mismatch inside it
            at = r.randrange(len(w))
            w = w[:at] + r.choice("ACGTN") + w[at + 1:]
        return w

    records = []
    for j in range(r.randint(1, 8)):
        kind = r.random()
        if kind < 0.12:
            seq = r.choice(list(configured.values()))           # an identical copy
            if r.random() < 0.3:
                seq = seq.lower()
        elif kind < 0.22:
            seq = _rand(r, r.randint(0, k - 1))                 # shorter than the window
        else:
            parts = [_rand(r, r.randint(0, 60))]
            for _ in range(r.randint(0, 4)):
                parts += [planted(), _rand(r, r.randint(0, 3 * k), r.choice(["ACGT", "AC", "A",
                                                                             "ACGTN"]))]
            seq = "".join(parts)
        tid = "ENST%011d" % r.choice([1, 2, 3, 500 + j, 500 + j, 900 + j])
        head = r.choice(["%s.%d|ENSG01|-|-|X-201|X|%d|protein_coding|" % (tid, j + 1, len(seq)),
                         "%s cdna chromosome:GRCh38:1:1:10:1 gene:ENSG01" % tid, tid])
        width = r.choice([60, 61, 7, 1000])
        lines = [seq[a:a + width] for a in range(0, len(seq), width)] or ([""] if r.random()
                                                                         < 0.5 else [])
        if r.random() < 0.2:
            lines = [ln + r.choice([" ", "\r", "  \t"]) for ln in lines]
        records.append(">%s\n%s" % (head, "".join(ln + "\n" for ln in lines)))
    text = "".join(records)
    gz = r.random() < 0.3
    path = tmp_path / ("bg%d.fa%s" % (i, ".gz" if gz else ""))
    if gz:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    exclude = r.sample(sorted(configured), r.randint(0, len(configured)))
    identical = None
    if r.random() < 0.7:
        identical = {}
        for t in sorted(configured):
            identical.setdefault(configured[t].strip().upper(), t)
    return {"path": path, "query": query, "k": k, "canonical": canonical,
            "exclude_ids": exclude, "identical": identical, "seed": seed}


def _both(c):
    want_out, got_out = {}, {}
    want = _scan_241(c["path"], c["query"], c["k"], canonical=c["canonical"],
                     exclude_ids=c["exclude_ids"], identical=c["identical"],
                     identical_out=want_out)
    got = I.scan_background_fasta(c["path"], c["query"], c["k"], canonical=c["canonical"],
                                  exclude_ids=c["exclude_ids"], identical=c["identical"],
                                  identical_out=got_out, seed=c["seed"])
    return (want, want_out), (got, got_out)


def test_the_scan_gives_2_4_1s_answer_on_random_cases(tmp_path):
    r = random.Random(20261002)
    hits = 0
    for i in range(600):
        c = _case(r, tmp_path, i)
        want, got = _both(c)
        assert got == want, "case %d: %r" % (i, {key: c[key] for key in ("k", "canonical",
                                                                          "seed")})
        hits += bool(want[0])
    # the cases are not vacuous: most find something
    assert hits > 400


@pytest.mark.parametrize("seed", [1, 5, 16, 31, 40])
def test_a_window_in_either_orientation_at_every_offset_is_found(tmp_path, seed):
    """Every start position modulo the stride, so a stride one too long misses some."""
    r = random.Random(7)
    k = 31
    conf = _rand(r, 200)
    query = I.kmers(conf, k)
    path = tmp_path / "bg.fa"
    recs = []
    for off in range(40):
        w = conf[off:off + k]
        recs.append(">R%d\n%s%s%s\n" % (off, "A" * off,
                                        I.revcomp(w) if off % 2 else w, "C" * 7))
    path.write_text("".join(recs))
    want = _scan_241(path, query, k)
    assert len(want) == 40
    assert I.scan_background_fasta(path, query, k, seed=seed) == want


def _records_ref(path, query_kmers, k, canonical=True, exclude_ids=(), identical=None,
                 identical_out=None):
    """Per record, what 2.4.1's scan looked at: every window folded and looked up."""
    query = set(query_kmers)
    if not query:
        return []
    exclude = {str(i).split(".")[0] for i in exclude_ids}
    identical = identical or {}
    out = []

    def _consume(chunks, tid):
        if not chunks or tid in exclude:
            return
        seq = "".join(chunks).upper()
        if seq in identical:
            if identical_out is not None:
                identical_out[tid] = identical[seq]
            return
        hits = set()
        for i in range(len(seq) - k + 1):
            km = I.canonical_kmer(seq[i:i + k]) if canonical else seq[i:i + k]
            if km in query:
                hits.add(km)
        out.append((tid, seq, hits))

    chunks, tid = [], None
    with io.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                _consume(chunks, tid)
                head = line[1:].strip()
                tid = head.split("|")[0].split()[0].split(".")[0] if head else None
                chunks = []
            else:
                chunks.append(line.strip())
        _consume(chunks, tid)
    return out


def test_each_record_holds_what_the_simple_scan_finds_in_it(tmp_path):
    r = random.Random(20261003)
    for i in range(400):
        c = _case(r, tmp_path, i)
        args = (c["path"], c["query"], c["k"], c["canonical"], c["exclude_ids"], c["identical"])
        want_out, got_out = {}, {}
        want = _records_ref(*args, identical_out=want_out)
        got = list(I._scan_records(*args, identical_out=got_out, seed=c["seed"]))
        assert got == want and got_out == want_out, "case %d" % i


def test_the_competitors_are_the_records_that_share_a_window_and_the_ones_asked_for(
        tmp_path):
    r = random.Random(3)
    conf = _rand(r, 300)
    query = I.kmers(conf, 31)
    near, far, other = conf[40:140], _rand(r, 200), _rand(r, 150)
    path = tmp_path / "bg.fa"
    path.write_text(">A.1|g\n%s\n>B.2|g\n%s\n>C.1\n%s\n>A.4\n%s\n>A.5\n%s\n>\n%s\n"
                    % (near, far, far, near, other, I.revcomp(near)))
    got = I.scan_fasta_competitors(path, query, 31, keep_ids=["C.7"])
    assert list(got) == ["A", "C", "record6"]          # file order; the copy of A dropped
    assert got["A"] == (near, I.kmers(near, 31))
    assert got["C"] == (far, set())                    # kept because it was asked for
    assert got["record6"][1] == I.kmers(near, 31)
    # a second A with other sequence that shares a window is kept under its own name
    path.write_text(">A\n%s\n>A\n%s\n" % (near, conf[100:200]))
    assert list(I.scan_fasta_competitors(path, query, 31)) == ["A", "A#2"]


def test_an_empty_query_reads_nothing_as_2_4_1_did(tmp_path):
    path = tmp_path / "bg.fa"
    path.write_text(">A\nACGTACGTAC\n")
    want_out, got_out = {}, {}
    assert I.scan_background_fasta(path, set(), 5, identical={"ACGTACGTAC": "T"},
                                   identical_out=got_out) == set() \
        == _scan_241(path, set(), 5, identical={"ACGTACGTAC": "T"}, identical_out=want_out)
    assert got_out == want_out == {}
    assert I.scan_fasta_competitors(path, set(), 5, keep_ids=["A"]) == {}


def test_a_third_record_with_one_id_and_other_sequence_is_numbered_on(tmp_path):
    r = random.Random(4)
    conf = _rand(r, 300)
    query = I.kmers(conf, 31)
    path = tmp_path / "bg.fa"
    path.write_text(">A\n%s\n>A\n%s\n>A\n%s\n" % (conf[:100], conf[100:200], conf[200:]))
    assert list(I.scan_fasta_competitors(path, query, 31)) == ["A", "A#2", "A#3"]
