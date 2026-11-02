"""Several configs' passes over one background FASTA give each config 2.10.0's answer.

``identifiability._whole_index_many`` reads the FASTA once for the first pass of several
whole-index calls, and once more for the second pass; ``index_scope._fasta_copies_many``
reads the headers once for several same-name-copy calls.  Each call's answer must be the
one 2.10.0 gives it alone.  2.10.0's functions are in ``oracle_2_10``, copied verbatim; the
first test holds them to 2.10.0's text.

The cases are drawn at random from a fixed seed, and hold what the records of a FASTA can
be to a call: a copy of a configured or a gene-background sequence, of an earlier record,
or of a sequence shorter than the window; a gene-background transcript, with the gene's
sequence or another; a decoy; a record longer than ``LONG_RECORD``; a repeat many records
hold; and one config's transcript as another's gene background or outside record.
``test_the_cases_hold_every_kind`` checks that every kind is drawn.
"""
import hashlib
import inspect
import random
from collections import Counter

import oracle_2_10 as O
import pytest

from isoform_dominance import identifiability as I
from isoform_dominance import index_scope

#: SHA-256 of each function's text at commit 4e082e1 (v2.10.0), as ``git show`` gives it.
ORACLE_SHA256 = {
    "_scan_records": "2352c3248cbd6806ba8342f7151c7c5dee325ceb4267eeb0efc8f965b31510a2",
    "_Names": "a08d3e67044f61de9fa4e2a42505d55425e0546f735d1bc06585d76e7a552280",
    "_whole_index": "b3e93368116c2f7a4b039a36b1c11be3b5c0b9dbe502cf0097f515c4c4a63e0b",
    "_relevant": "5d0e0bd166d7da9a0c504b01a9e53f52f4afda62ce2ed32be60cbfe618e44984",
    "copies_in": "82a2c7730c5377c28dd52fbc030d87ed1af19041920f9c68fab30d50408c6dc9",
    "fasta_copies": "de6554065029786d3800ad62740dcd58f63e320740b5f33836e4c1dd4952b862",
}

N_WHOLE_CASES = 320
N_COPIES_CASES = 320
#: the package's LONG_RECORD, before a test lowers it: a case's own does not hang on the last
LONG_RECORD = I.LONG_RECORD


def test_the_oracle_is_2_10_0s_text():
    got = {n: hashlib.sha256(inspect.getsource(getattr(O, n)).encode()).hexdigest()
           for n in ORACLE_SHA256}
    assert got == ORACLE_SHA256


# --------------------------------------------------------------------------- #
# the cases
# --------------------------------------------------------------------------- #
def _rand(r, n, alphabet="ACGT"):
    return "".join(r.choice(alphabet) for _ in range(n))


def _mutate(r, s, n):
    s = list(s)
    for _ in range(min(n, len(s))):
        i = r.randrange(len(s))
        s[i] = r.choice("ACGT".replace(s[i], "") or "A")
    return "".join(s)


def _tid(n):
    return "ENST%011d" % n


def _header(r, tid, gene_id, name, kinds):
    """A header for ``tid``: GENCODE's, Ensembl's or the id alone, with a version or not."""
    shown = tid + (".%d" % r.randint(1, 9) if r.random() < 0.7 else "")
    if shown != tid:
        kinds["versioned_id"] += 1
    form = r.choice(["gencode", "ensembl", "bare"])
    kinds["header_" + form] += 1
    if form == "gencode":
        return "%s|%s.1|-|-|%s-201|%s|99|protein_coding|" % (shown, gene_id, name, name)
    if form == "ensembl":
        return "%s cdna chromosome:GRCh38:1:1:99:1 gene:%s.4 gene_symbol:%s" % (
            shown, gene_id, name)
    return shown


def _lines(r, seq, kinds):
    """``seq`` as FASTA lines: wrapped or not, some with white space around them."""
    if not seq:
        return [""]
    width = r.choice([len(seq), 7, 13, 60])
    out = [seq[i:i + width] for i in range(0, len(seq), width)]
    if r.random() < 0.15:
        out = [" " + x + " \t" for x in out]
        kinds["white_space_in_lines"] += 1
    return out


def _column_args(r, window, canonical, needed, seqs, gene, keep_duplicates,
                 max_window_records, kinds):
    """The arguments ``analyze`` gives ``_whole_index`` for one config (less the path)."""
    configured = {}
    for t in sorted(set(needed)):
        configured.setdefault(seqs[t].strip().upper(), t)
    cols, first = {}, {}
    for t in sorted(gene):
        n = gene[t].strip().upper()
        if not keep_duplicates and n in configured:
            continue
        if not keep_duplicates and n in first:
            continue
        first[n] = t
        cols[t] = gene[t]
    configured_windows = set()
    for t in needed:
        configured_windows.update(I.kmer_track(seqs[t], window, canonical))
    column_windows = set(configured_windows)
    for s in cols.values():
        column_windows.update(I.kmer_track(s, window, canonical))
    if column_windows and r.random() < 0.15:
        # windows _scan_records never finds: the wrong length, and not canonical
        w = sorted(column_windows)[0]
        column_windows |= {w + "A", w[:-1]}
        if canonical:
            column_windows.add(I.revcomp(w) if I.revcomp(w) != w else w)
        kinds["odd_query_windows"] += 1
    if not column_windows:
        kinds["empty_column_windows"] += 1
    return {"column_windows": column_windows, "configured_windows": configured_windows,
            "gene_ids": set(gene), "gene_sequences": first, "configured": configured,
            "exclude_ids": list(needed), "keep_duplicates": keep_duplicates,
            "max_window_records": max_window_records}


def _whole_case(seed, tmp_path):
    """``(path, window, canonical, decoys, long_record, [per-config arguments], kinds)``."""
    r = random.Random(seed)
    kinds = Counter()
    window = r.choice([4, 6, 9, 13, 15, 16, 17, 21])
    canonical = r.random() < 0.75
    if window < I.FASTA_SEED:
        kinds["window_below_seed"] += 1
    if not canonical:
        kinds["canonical_false"] += 1
    exons = [_rand(r, r.randint(window + 3, 3 * window + 20)) for _ in range(9)]
    rep = _rand(r, window + r.randint(1, 2 * window))
    short = [_rand(r, r.randint(1, window - 1)) for _ in range(2)]

    def transcript():
        s = "".join(r.sample(exons, r.randint(1, 3)))
        return s + rep if r.random() < 0.3 else s

    n_jobs = r.randint(2, 5)
    ids = iter(range(1, 10_000))
    jobs = []                                   # (needed, seqs, gene)
    for _ in range(n_jobs):
        seqs = {_tid(next(ids)): transcript() for _ in range(r.randint(1, 3))}
        gene = {_tid(next(ids)): transcript() for _ in range(r.randint(0, 3))}
        u = r.random()
        if u < 0.12:
            # every column shorter than the window: no column window
            seqs = {t: r.choice(short) for t in seqs}
            gene = {t: r.choice(short) for t in gene}
        elif u < 0.25:
            seqs[next(iter(seqs))] = r.choice(short)
            kinds["short_configured"] += 1
            if gene:
                gene[next(iter(gene))] = r.choice(short)
                kinds["short_gene_column"] += 1
        elif u < 0.3:
            seqs[next(iter(seqs))] = ""
            kinds["empty_configured"] += 1
        jobs.append([list(seqs), seqs, gene])
    if r.random() < 0.5:
        # one config's transcript in another's gene background
        a, b = r.sample(range(n_jobs), 2)
        t = jobs[a][0][0]
        jobs[b][2][t] = jobs[a][1][t]
        kinds["configured_is_anothers_gene"] += 1

    records = []                                # (header or None, lines)
    gene_name = {i: "GENE%d" % i for i in range(n_jobs)}

    def add(tid, seq, kind, job=None):
        name = gene_name[job] if job is not None else r.choice(["OTHER", "GENE0"])
        gid = "ENSG%011d" % (job if job is not None else 900 + r.randint(0, 3))
        records.append((_header(r, tid, gid, name, kinds), _lines(r, seq, kinds)))
        kinds[kind] += 1

    if r.random() < 0.15:
        records.append((None, _lines(r, transcript(), kinds)))
        kinds["lines_before_first_header"] += 1
    decoys = []
    for _ in range(r.randint(4, 24)):
        j = r.randrange(n_jobs)
        needed, seqs, gene = jobs[j]
        kind = r.choices(
            ["paralogue", "repeat", "second", "copy_configured", "copy_gene", "copy_earlier",
             "short_copy", "gene_id_other", "gene_id_same", "id_twice", "no_id", "blank",
             "no_lines", "decoy", "lower_n_iupac", "configured_id", "random"],
            [4, 4, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2])[0]
        tid = _tid(next(ids))
        if kind == "paralogue":
            add(tid, _mutate(r, r.choice(exons), r.randint(0, 2)) + _rand(r, 20), kind, j)
        elif kind == "repeat":
            for _ in range(r.randint(2, 6)):
                add(_tid(next(ids)), _rand(r, 10) + rep + _rand(r, 12), kind)
        elif kind == "second" and records and "".join(records[-1][1]).strip():
            # shares the start of the record before: an added record's own window
            add(tid, "".join(x.strip() for x in records[-1][1])[:window + 4] + _rand(r, 25),
                kind)
        elif kind == "copy_configured":
            add(tid, seqs[r.choice(needed)], kind, j)
        elif kind == "copy_gene" and gene:
            add(tid, gene[r.choice(sorted(gene))], kind, j)
        elif kind == "copy_earlier" and records and records[-1][0]:
            add(tid, "".join(x.strip() for x in records[-1][1]), kind)
        elif kind == "short_copy":
            add(tid, r.choice(short), kind, j)
        elif kind == "gene_id_other" and gene:
            add(r.choice(sorted(gene)), transcript(), kind, j)
        elif kind == "gene_id_same" and gene:
            t = r.choice(sorted(gene))
            add(t, gene[t], kind, j)
        elif kind == "id_twice" and records and (records[-1][0] or "").split("|")[0].split():
            first = records[-1][0].split("|")[0].split()[0].split(".")[0]
            add(first, transcript() if r.random() < 0.5 else
                "".join(x.strip() for x in records[-1][1]), kind)
        elif kind == "no_id":
            records.append((r.choice([" ", "|x", "", "|ENSG00000000001.1|-"]),
                            _lines(r, transcript(), kinds)))
            kinds[kind] += 1
        elif kind == "blank":
            records.append((tid, ["", "  "]))
            kinds[kind] += 1
        elif kind == "no_lines":
            records.append((tid, []))
            kinds[kind] += 1
        elif kind == "decoy":
            name = "chrD%d" % len(decoys)
            decoys.append(name)
            by_field = r.random() < 0.5
            records.append(((name + "|extra words") if by_field else (name + " more"),
                            _lines(r, transcript() + rep, kinds)))
            kinds["decoy_by_field" if by_field else "decoy_by_word"] += 1
        elif kind == "lower_n_iupac":
            s = list(transcript())
            for _ in range(r.randint(1, 4)):
                s[r.randrange(len(s))] = r.choice("NRYKMnacgt")
            add(tid, "".join(s), kind, j)
        elif kind == "configured_id":
            # one config's transcript is a record: the others meet it as an outside record
            add(needed[0], seqs[needed[0]] if r.random() < 0.5 else transcript(), kind)
        else:
            add(tid, transcript(), "random")
    if decoys and r.random() < 0.3:
        decoys.append("chrAbsent")
    long_record = LONG_RECORD
    lengths = sorted(len("".join(x.strip() for x in lines)) for _, lines in records)
    if r.random() < 0.3 and lengths[-1] > 1:
        long_record = lengths[-1] - 1 if r.random() < 0.5 else lengths[len(lengths) // 2]
        kinds["long_record"] += 1

    text = []
    for head, lines in records:
        if head is not None:
            text.append(">" + head + "\n")
        text.extend(x + "\n" for x in lines)
    path = tmp_path / ("whole%d.fa" % seed)
    path.write_text("".join(text))

    args = []
    for needed, seqs, gene in jobs:
        keep = r.random() < 0.25
        if keep:
            kinds["keep_duplicates"] += 1
        m = r.choice([0, 1, 20])
        kinds["max_window_records_%d" % m] += 1
        args.append(_column_args(r, window, canonical, needed, seqs,
                                 {t: s for t, s in gene.items() if t not in needed},
                                 keep, m, kinds))
    return str(path), window, canonical, decoys, long_record, args, kinds


WHOLE_KINDS = [
    "window_below_seed", "canonical_false", "short_configured", "short_gene_column",
    "empty_configured", "configured_is_anothers_gene", "lines_before_first_header",
    "paralogue", "repeat", "second", "copy_configured", "copy_gene", "copy_earlier",
    "short_copy", "gene_id_other", "gene_id_same", "id_twice", "no_id", "blank", "no_lines",
    "decoy_by_field", "decoy_by_word", "lower_n_iupac", "configured_id", "long_record",
    "keep_duplicates", "max_window_records_0", "max_window_records_1",
    "max_window_records_20", "empty_column_windows", "odd_query_windows", "versioned_id",
    "header_gencode", "header_ensembl", "header_bare", "white_space_in_lines",
]


def _whole_cases(tmp_path):
    for seed in range(N_WHOLE_CASES):
        yield _whole_case(seed, tmp_path)


def _long(monkeypatch, n):
    monkeypatch.setattr(I, "LONG_RECORD", n)
    monkeypatch.setattr(O, "LONG_RECORD", n)


# --------------------------------------------------------------------------- #
# whole index
# --------------------------------------------------------------------------- #
def test_the_cases_hold_every_kind(tmp_path):
    total = Counter()
    for *_, kinds in _whole_cases(tmp_path):
        total.update(kinds)
    assert {k: total[k] for k in WHOLE_KINDS if total[k] < 10} == {}


def test_many_calls_read_together_give_each_its_own_answer(tmp_path, monkeypatch):
    checked = Counter()
    for path, window, canonical, decoys, long_record, args, _ in _whole_cases(tmp_path):
        _long(monkeypatch, long_record)
        got = I._whole_index_many(path, args, window, canonical, decoys)
        assert len(got) == len(args)
        for kw, mine in zip(args, got, strict=True):
            want = O._whole_index(path, window=window, canonical=canonical, decoys=decoys,
                                  **kw)
            assert mine == want
            checked["calls"] += 1
            checked["added"] += len(want["added"])
            checked["left_out"] += want["left_out"]
            checked["second_round"] += bool(want["gene_records"])
            checked["read"] += want["stats"]["read"]
    # the cases reach every branch that matters
    assert min(checked.values()) > 50, checked


def test_the_packages_whole_index_is_2_10_0s(tmp_path, monkeypatch):
    for path, window, canonical, decoys, long_record, args, _ in _whole_cases(tmp_path):
        _long(monkeypatch, long_record)
        for kw in args:
            assert I._whole_index(path, window=window, canonical=canonical, decoys=decoys,
                                  **kw) == O._whole_index(path, window=window,
                                                          canonical=canonical, decoys=decoys,
                                                          **kw)


def test_one_call_read_with_others_or_alone_is_the_same(tmp_path, monkeypatch):
    path, window, canonical, decoys, long_record, args, _ = _whole_case(7, tmp_path)
    _long(monkeypatch, long_record)
    together = I._whole_index_many(path, args, window, canonical, decoys)
    alone = [I._whole_index_many(path, [kw], window, canonical, decoys)[0] for kw in args]
    assert together == alone
    assert I._whole_index_many(path, [], window, canonical, decoys) == []


def test_a_call_that_cannot_be_handled_is_that_calls_error(tmp_path, monkeypatch):
    path, window, canonical, decoys, long_record, args, _ = _whole_case(3, tmp_path)
    _long(monkeypatch, long_record)
    bad = dict(args[0], max_window_records=None)    # m(w) <= None: a TypeError
    bad["configured_windows"] = set(bad["column_windows"])
    got = I._whole_index_many(path, [bad] + args, window, canonical, decoys)
    with pytest.raises(TypeError) as want:
        O._whole_index(path, window=window, canonical=canonical, decoys=decoys, **bad)
    assert type(got[0]) is TypeError and got[0].args == want.value.args
    assert got[1:] == [O._whole_index(path, window=window, canonical=canonical,
                                      decoys=decoys, **kw) for kw in args]


# --------------------------------------------------------------------------- #
# same-name copies
# --------------------------------------------------------------------------- #
def _copies_case(seed, tmp_path):
    """``(path, [(target ids, gene names, versions or None)], kinds)``."""
    r = random.Random(1000 + seed)
    kinds = Counter()
    genes = []                                   # (gene id, name)
    for i in range(r.randint(2, 6)):
        genes.append(("ENSG%011d" % i, r.choice(["ALPHA", "BETA", "GAMMA", "DELTA"])))
    heads, ids = [], iter(range(1, 10_000))
    by_gene = {}
    for _ in range(r.randint(3, 30)):
        gid, name = r.choice(genes)
        t = _tid(next(ids))
        by_gene.setdefault(gid, []).append(t)
        shown = t + (".%d" % r.randint(1, 5) if r.random() < 0.7 else "")
        form = r.choices(["gencode", "ensembl_ref", "ensembl_off", "bare", "broken"],
                         [4, 2, 2, 1, 1])[0]
        kinds["header_" + form] += 1
        if form == "gencode":
            heads.append("%s|%s.2|-|-|%s-201|%s|99|protein_coding|" % (shown, gid, name, name))
        elif form.startswith("ensembl"):
            region = ("chromosome:GRCh38:%s:1:9:1" % r.choice(["1", "X", "Y"])
                      if form == "ensembl_ref" else "scaffold:GRCh38:KI270713.1:1:9:1")
            heads.append("%s cdna %s gene:%s.3 gene_symbol:%s" % (shown, region, gid, name))
        elif form == "bare":
            heads.append(shown)
        else:
            heads.append(r.choice(["", " ", "|x", "%s|%s||" % (shown, gid)]))
        if r.random() < 0.1:
            heads.append(heads[-1])
            kinds["header_twice"] += 1
    requests = []
    for _ in range(r.randint(1, 5)):
        gid, name = r.choice(genes)
        mine = by_gene.get(gid, [])
        targets = r.sample(mine, min(len(mine), r.randint(1, 3))) if mine else []
        if r.random() < 0.3:
            targets.append(_tid(9_000 + r.randint(0, 5)))         # not in the FASTA
            kinds["target_absent"] += 1
        targets = [t + ".7" if r.random() < 0.3 else t for t in targets]
        u = r.random()
        if u < 0.35:
            names = [name]
            kinds["names_the_genes"] += 1
        elif u < 0.6:
            names = [r.choice(["ALPHA", "BETA", "GAMMA", "DELTA", "OTHER"])]
            kinds["names_maybe_another"] += 1
        elif u < 0.8:
            names = ()
            kinds["names_none"] += 1
        elif u < 0.9:
            names = ["", None]
            kinds["names_empty"] += 1
        else:
            names = [["a list"]]                    # a config whose gene is a list
            kinds["names_unhashable"] += 1
        requests.append((targets, names, {} if r.random() < 0.8 else None))
    path = tmp_path / ("copies%d.fa" % seed)
    path.write_text("".join(">%s\nACGT\n" % h for h in heads))
    return str(path), requests, kinds


COPIES_KINDS = ["header_gencode", "header_ensembl_ref", "header_ensembl_off", "header_bare",
                "header_broken", "header_twice", "target_absent", "names_the_genes",
                "names_maybe_another", "names_none", "names_empty", "names_unhashable"]


def _outcome(call, *args, **kwargs):
    try:
        return "ok", call(*args, **kwargs)
    except Exception as e:
        return type(e), e.args


def test_many_copy_calls_from_one_read_of_the_headers_give_each_its_answer(tmp_path):
    total, checked = Counter(), Counter()
    for seed in range(N_COPIES_CASES):
        path, requests, kinds = _copies_case(seed, tmp_path)
        total.update(kinds)
        mine_versions = [None if v is None else {} for _, _, v in requests]
        got = index_scope._fasta_copies_many(
            path, [(t, n, mv) for (t, n, _), mv in zip(requests, mine_versions, strict=True)])
        for (targets, names, versions), mine, mv in zip(requests, got, mine_versions,
                                                        strict=True):
            want = _outcome(O.fasta_copies, path, targets, names, versions=versions)
            if want[0] == "ok":
                assert mine == want[1]
                checked["copies"] += bool(want[1])
            else:
                assert (type(mine), mine.args) == want
                checked["errors"] += 1
            assert mv == versions
            if versions is not None:
                assert list(mv) == list(versions)           # in the same order
                checked["versions"] += bool(versions)
    assert {k: total[k] for k in COPIES_KINDS if total[k] < 10} == {}
    assert min(checked.values()) > 20, checked


def test_copy_calls_whose_name_differs_need_the_headers_again_and_get_them(tmp_path):
    # the copy comes first, and the config names its gene otherwise: copies_in reads twice
    path = tmp_path / "two.fa"
    path.write_text(">ENST00000000009.1|ENSG00000000009.1|-|-|X-201|SMN|9|pc|\nACGT\n"
                    ">ENST00000000001.2|ENSG00000000001.1|-|-|X-201|SMN|9|pc|\nACGT\n")
    requests = [(["ENST00000000001"], ["SMN1"], {}), (["ENST00000000001"], (), {})]
    got = index_scope._fasta_copies_many(str(path), requests)
    want = [O.fasta_copies(str(path), t, n, versions={}) for t, n, _ in requests]
    assert got == want
    assert got[0] == [{"gene_id": "ENSG00000000009", "gene_name": "SMN", "region": None,
                       "transcripts": ["ENST00000000009"]}]
