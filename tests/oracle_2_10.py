"""The functions 2.11.0's shared FASTA reads are held to: v2.10.0's, verbatim.

Copied character for character from commit 4e082e1 (the 2.10.0 release, tag v2.10.0):

- ``src/isoform_dominance/identifiability.py``: ``_scan_records``, ``_Names``,
  ``_whole_index``;
- ``src/isoform_dominance/index_scope.py``: ``_relevant``, ``copies_in``,
  ``fasta_copies``.

What they call is imported from the package, which this round does not change:
``kmer_track``, ``canonical_kmer``, ``revcomp``, ``index_scope.parse_header``,
``index_scope.same_name_copies``, ``index_scope._unversioned``, ``io.open_text`` and the
constants.  ``LONG_RECORD`` is this module's own, read when called, as the package reads
its own: a test that lowers it lowers both.  ``tests/test_whole_index_many.py`` holds each
function's text to the SHA-256 of 2.10.0's.  The name does not start with ``test_``, so
pytest does not collect it.
"""
import hashlib
from collections import Counter

from isoform_dominance import io
from isoform_dominance.identifiability import (DEFAULT_MAX_WINDOW_RECORDS, FASTA_SEED,
                                               LONG_RECORD, canonical_kmer, kmer_track,
                                               revcomp)
from isoform_dominance.index_scope import _unversioned, parse_header, same_name_copies
from isoform_dominance.io import open_text


def _scan_records(path, query_kmers, k, canonical=True, exclude_ids=(), identical=None,
                  identical_out=None, seed=FASTA_SEED, decoys=(), stats=None, every=False):
    """Yield ``(record id, sequence, windows)`` for each record of a background FASTA.

    The sequence is the record's lines, stripped, joined and upper-cased; ``windows`` is
    the set of ``query_kmers`` the record holds, in the form ``query_kmers`` gives them.
    Records that :func:`scan_background_fasta` skips -- excluded, identical to a
    configured transcript, or with no sequence -- are not yielded, and nothing is
    yielded for an empty query.  See :func:`scan_background_fasta` for the arguments and
    for how a window is found.  The id is None when the header gives none.

    With ``every``, a record identical to a configured transcript is yielded too, and each
    record as ``(place, id, sequence, windows, same)``: ``place`` counts the headers of
    the file from 1, ``same`` is the configured transcript the record is identical to
    (``windows`` then None), or None.

    ``decoys`` are record names, as Salmon's ``decoys.txt`` gives them: the header's
    first whitespace-delimited word, or that word's first ``|``-delimited field.  Such a
    record is skipped unread.  ``stats``, a dict, receives ``decoys_skipped`` (records),
    ``decoys_found`` (the set of names met) and ``long_records``, ``[(id, length)]`` of
    the records read that are longer than :data:`LONG_RECORD`.
    """
    decoys = frozenset(decoys)
    if stats is not None:
        stats.update(read=False, decoys_skipped=0, decoys_found=set(), long_records=[])
    query = set(query_kmers)
    if not query:
        return
    if stats is not None:
        stats["read"] = True
    exclude = {str(i).split(".")[0] for i in exclude_ids}
    identical = identical or {}
    lengths = {len(s) for s in identical}

    # each way a record's window can read as a query window, to the query window it is
    look = {}
    for q in query:
        if len(q) != k:
            continue                    # a window of the record is k long
        if canonical:
            if canonical_kmer(q) != q:
                continue                # a folded window is canonical, so never this
            look[revcomp(q)] = q
        look[q] = q
    s = max(1, min(k, seed))
    t = max(1, k - s + 1)
    seeds = {w[j:j + s] for w in look for j in range(t)}

    def _hits(seq):
        hits = set()
        n = len(seq)
        if k <= 0 or n < k:
            return hits
        for p in range(0, n - s + 1, t):
            if seq[p:p + s] in seeds:
                # the windows whose seed position is p: they start at p - t + 1 .. p
                for i in range(max(0, p - t + 1), min(p, n - k) + 1):
                    hit = look.get(seq[i:i + k])
                    if hit is not None:
                        hits.add(hit)
        return hits

    def _consume(chunks, tid, place):
        if not chunks or tid in exclude:
            return None
        seq = "".join(chunks).upper()
        if len(seq) in lengths and seq in identical:
            if identical_out is not None:
                identical_out[tid] = identical[seq]
            return (place, tid, seq, None, identical[seq]) if every else None
        if not seq:
            return None             # blank lines only: no sequence, so no window either
        if stats is not None and len(seq) > LONG_RECORD:
            stats["long_records"].append((tid, len(seq)))
        return (place, tid, seq, _hits(seq), None) if every else (tid, seq, _hits(seq))

    chunks, tid, skip, place = [], None, False, 0
    with io.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                rec = None if skip else _consume(chunks, tid, place)
                if rec is not None:
                    yield rec
                place += 1
                head = line[1:].strip()
                # the first word of the first |-field, without its version; a header that
                # starts with "|" has none (2.4.1 raised IndexError)
                first = head.split("|")[0].split()
                tid = first[0].split(".")[0] if first else None
                chunks = []
                if decoys:
                    word = head.split()[0] if head else ""
                    name = next((x for x in (word, word.split("|")[0]) if x in decoys), None)
                    skip = name is not None
                    if skip and stats is not None:
                        stats["decoys_skipped"] += 1
                        stats["decoys_found"].add(name)
            elif not skip:
                chunks.append(line.strip())
        rec = None if skip else _consume(chunks, tid, place)
        if rec is not None:
            yield rec


class _Names:
    """The name each reported record goes by; see :func:`scan_fasta_competitors`.

    ``names(id, place, sequence)`` gives the name, or None for a record with an id and a
    sequence already reported together.  Sequences are kept as digests.
    """

    def __init__(self):
        self.taken = set()
        self.seen = {}                  # id -> digests of the sequences reported under it

    def __call__(self, tid, place, seq):
        rid = tid if tid else "record%d" % place
        digest = hashlib.blake2b(seq.encode(), digest_size=16).digest()
        mine = self.seen.setdefault(rid, set())
        if digest in mine:
            return None
        mine.add(digest)
        name, m = rid, 2
        while name in self.taken:
            name, m = "%s#%d" % (rid, m), m + 1
        self.taken.add(name)
        return name

    def claim(self, tid, seq):
        """Reserve ``tid`` for a column with this sequence."""
        self.taken.add(tid)
        self.seen.setdefault(tid, set()).add(
            hashlib.blake2b(seq.encode(), digest_size=16).digest())


def _whole_index(path, column_windows, configured_windows, gene_ids, gene_sequences,
                 configured, window, canonical=True, exclude_ids=(), keep_duplicates=False,
                 max_window_records=DEFAULT_MAX_WINDOW_RECORDS, decoys=(),
                 seed=FASTA_SEED):
    """The records of a background FASTA that join the system, and the windows to drop.

    ``column_windows`` are the windows of the columns (the configured transcripts and the
    gene background), ``configured_windows`` those of the configured transcripts alone.
    The first record with a ``gene_ids`` id is that transcript of the gene background, and
    its sequence is returned, not judged.  Unless ``keep_duplicates``, a record with a
    configured transcript's sequence (``configured``, upper-cased sequence -> id), a gene
    column's (``gene_sequences``, the same) or an earlier record's sequence is a copy, not
    a record.  Every other record, not excluded and not a decoy, is an *outside* record.

    Pass 1 counts, for each column window ``w``, the outside records that hold it,
    ``m(w)``, each record once.  An outside record that holds a configured window with
    ``m(w) <= max_window_records`` is *added*, unless it is longer than
    :data:`LONG_RECORD`; every other outside record is *left out*.  Pass 2 looks for the
    windows of the added records in the left-out records.  Every window of a column or
    of an added record that a left-out record holds is *dropped*.

    Returns a dict: ``gene_records`` (id -> sequence), ``identical_to_configured`` and
    ``identical_to_background`` (record name -> the id it equals), ``added`` (name ->
    ``(sequence, number of distinct configured windows it holds)``), ``dropped`` (a set
    of windows), ``left_out`` (how many outside records that hold a window of the system
    were left out) and ``stats`` (see :func:`_scan_records`).  Record names are as
    :func:`scan_fasta_competitors` gives them.
    """
    stats = {}
    exclude = list(exclude_ids)
    _name = _Names()

    gene_records, same_conf, same_bg = {}, {}, {}
    gene_lengths = {len(x) for x in gene_sequences}
    first_of = {}                       # digest of a sharing record's sequence -> its name
    sharing = []                        # (name, place, hits, configured hits, sequence)
    m = Counter()
    settled = set()                     # places that are no outside record: columns, copies
    for place, tid, seq, hits, same in _scan_records(
            path, column_windows, window, canonical, exclude, None if keep_duplicates
            else configured, None, seed, decoys, stats, every=True):
        if tid in gene_ids and tid not in gene_records:
            gene_records[tid] = seq
            _name.claim(tid, seq)       # the column's name: a later record of it is <id>#2
            settled.add(place)
            continue
        if same is not None:
            rid = _name(tid, place, seq)
            if rid is not None:
                same_conf[rid] = same
            settled.add(place)
            continue
        if not keep_duplicates and len(seq) in gene_lengths and seq in gene_sequences:
            rid = _name(tid, place, seq)
            if rid is not None:
                same_bg[rid] = gene_sequences[seq]
            settled.add(place)
            continue
        if not hits:
            continue                    # outside, and holds no column window
        digest = None
        if not keep_duplicates:
            digest = hashlib.blake2b(seq.encode(), digest_size=16).digest()
            if digest in first_of:
                rid = _name(tid, place, seq)
                if rid is not None:
                    same_bg[rid] = first_of[digest]
                settled.add(place)
                continue
        rid = _name(tid, place, seq)
        if rid is None:
            settled.add(place)
            continue
        if digest is not None:
            first_of[digest] = rid
        mine = hits & configured_windows
        for w in hits:
            m[w] += 1
        # m only grows, so a record none of whose configured windows is within the limit
        # now never will be: its sequence need not be kept
        hopeful = len(seq) <= LONG_RECORD and any(m[w] <= max_window_records for w in mine)
        sharing.append((rid, place, hits, mine, seq if hopeful else None, digest))

    added = {}
    within = Counter()                  # per window, the added records that hold it
    left_places = set()
    left_seqs = set()                   # digests of the left-out records, copies counted once
    for rid, place, hits, mine, seq, digest in sharing:
        if seq is not None and any(m[w] <= max_window_records for w in mine):
            added[rid] = (seq, len(mine))
            within.update(hits)
            settled.add(place)
        else:
            left_places.add(place)
            left_seqs.add(digest if digest is not None else place)
    # a column window held by more outside records than were added is held by a left-out one
    dropped = {w for w, n in m.items() if n > within[w]}

    # pass 2: the added records' own windows, in the records left out
    query = set()
    for seq, _ in added.values():
        query.update(kmer_track(seq, window, canonical))
    query -= column_windows             # those were settled in pass 1
    if query:
        for place, _, seq, hits, _ in _scan_records(
                path, query, window, canonical, exclude, None, None, seed, decoys,
                None, every=True):
            if hits and place not in settled:
                dropped |= hits         # a pass-1 left-out record's too
                if place not in left_places:
                    left_places.add(place)
                    # a copy is no record: each sequence is one record left out
                    left_seqs.add(place if keep_duplicates else
                                  hashlib.blake2b(seq.encode(), digest_size=16).digest())
    # a copy of a record left out is no column's copy: it is not reported as one
    columns = set(gene_sequences.values()) | set(added)
    same_bg = {rid: t for rid, t in same_bg.items() if t in columns}
    return {"gene_records": gene_records, "identical_to_configured": same_conf,
            "identical_to_background": same_bg, "added": added, "dropped": dropped,
            "left_out": len(left_seqs), "stats": stats}


def _relevant(heads, targets, names, versions=None):
    """The parsed ``heads`` that are a target transcript or carry one of ``names``; and in
    ``versions``, each target's versioned ids, the header's first ``|``-field or word."""
    out = []
    for head in heads:
        if versions is not None:
            first = head.strip().split("|")[0].split()
            if first and _unversioned(first[0]) in targets:
                versions.setdefault(_unversioned(first[0]), set()).add(first[0])
        r = parse_header(head)
        if r and (r["transcript"] in targets or r["gene_name"] in names):
            out.append(r)
    return out


def copies_in(read_heads, target_ids, gene_names=(), versions=None):
    """:func:`same_name_copies` over the headers ``read_heads()`` yields, holding only the
    records that can matter rather than every header of a whole-transcriptome FASTA.

    ``read_heads`` is called again, once, when the name the headers give the target
    transcripts is not among ``gene_names``: a copy can come before its reference gene.
    ``versions``, a dict, receives the versioned ids the headers give each target
    transcript (``{id: {versioned id, ...}}``).
    """
    targets = {_unversioned(t) for t in target_ids}
    names = {n for n in gene_names or () if n}
    recs = _relevant(read_heads(), targets, names, versions)
    found = {r["gene_name"] for r in recs if r["transcript"] in targets}
    if found - names:
        names |= found
        recs = _relevant(read_heads(), targets, names)
    return same_name_copies(recs, targets, names)


def fasta_copies(path, target_ids, gene_names=(), versions=None):
    """:func:`copies_in` for a plain or gzipped FASTA (told by its magic bytes)."""
    def heads():
        with open_text(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    yield line[1:]
    return copies_in(heads, target_ids, gene_names, versions)
