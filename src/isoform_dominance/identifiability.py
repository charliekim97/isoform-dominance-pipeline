"""Group-level short-read identifiability for functional isoform classes.

Short-read quantifiers apportion fragments among transcripts by solving a linear
inverse problem: the expected count of each observable fragment class is a linear
function of the transcript abundances.  A quantity is *estimable* from that system
only when its coefficient vector lies in the row space of the design matrix, and it
is estimable *usefully* only when the corresponding conditioning factor is small.

Transcript-level identifiability of this system has been studied before
(Hiller et al. 2009, doi:10.1093/bioinformatics/btp544; Ferrer-Bonsoms et al. 2022,
doi:10.1093/bioinformatics/btab873), and ``terminus`` (Sarkar et al. 2020,
doi:10.1093/bioinformatics/btaa448) discovers, post hoc and from the data, transcript
groups whose *totals* carry low inferential uncertainty.

This module answers the question that sits in front of those: given two isoform
classes that the user has defined **on biological grounds, before quantification**,
is the class total -- and the contrast between the two class totals -- estimable from
short reads at all, and with what precision?  The relevant estimands are

    s_g = sum_{t in g} theta_t          (one per group)
    d   = s_A - s_B                     (the primary comparison)

and each is checked by the textbook estimability condition together with a structural
conditioning factor.  Three layers are reported, cheapest first:

1. **Sequence uniqueness** -- group-unique k-mers, the positions they cover, and how
   those positions are arranged into blocks.  Uniqueness is assessed against a
   *background* (by default every other transcript of the same gene, optionally an
   arbitrary FASTA such as the Salmon index), not merely against the other groups in
   the config.  Assessing uniqueness only among the configured groups overstates
   separability whenever a paralogue, a retained-intron transcript, or an unlisted
   isoform of the same gene carries the same sequence.
2. **Expected informative reads** -- how many read pairs, at a stated depth, class
   abundance, read length and fragment-length distribution, actually land such that a
   sequenced end covers a group-unique k-mer.  This converts a sequence property into
   an experimental-design quantity, and with it a counting-noise floor on the
   log2 class ratio.
3. **Estimability** -- rank, row-space residual and a structural conditioning factor
   for the class-collapsed compatibility system described above.  The conditioning
   factor is a geometry proxy, not a standard error; see :func:`estimability`.

All three are computed offline from sequence alone; nothing here needs the reads.
"""
import math
import urllib.parse
from collections import Counter
from urllib.error import HTTPError

import numpy as np

from . import ensembl, index_scope, io
from .ensembl import DEFAULT_RETRIES, DEFAULT_RETRY_WAIT

ENSEMBL = ensembl.SERVER

#: Default k-mer length; matches the Salmon index default.
DEFAULT_K = 31

#: Defaults for the read/fragment model used by :func:`expected_informative_reads`.
DEFAULT_READ_LENGTH = 100
DEFAULT_FRAG_MEAN = 200.0
DEFAULT_FRAG_SD = 60.0
DEFAULT_DEPTH = 30_000_000
DEFAULT_MEAN_EFFLEN = 1500.0
DEFAULT_TPM = 10.0

#: A group whose structural conditioning factor exceeds this is reported as weakly
#: identifiable rather than identifiable.  Provisional: chosen by reasoning about the
#: geometry, not calibrated against realised quantification error.  See
#: :func:`estimability`.
DEFAULT_CONDITIONING_TAU = 10.0

#: A class structurally estimable but expected to receive fewer than this many
#: unambiguously assignable fragments is downgraded to weakly identifiable: the
#: comparison is then limited by counting noise rather than by the design.
DEFAULT_MIN_INFORMATIVE_READS = 50.0

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


# --------------------------------------------------------------------------- #
# sequence retrieval
# --------------------------------------------------------------------------- #
def fetch_cdna(transcript_id, **retry):
    """Fetch one transcript's cDNA sequence from the Ensembl REST API.

    For more than one transcript use :func:`isoform_dominance.ensembl.fetch_cdna_batch`,
    which asks for 50 per request.  ``retry`` (``retries``, ``retry_wait``, ``timeout``)
    goes to :func:`isoform_dominance.ensembl.request`.
    """
    tid = transcript_id.split(".")[0]
    got = ensembl.fetch_cdna_batch([tid], **retry)
    if tid not in got:
        raise ValueError("Ensembl returned no cDNA for %s" % tid)
    return got[tid]


def fetch_gene(gene, species="homo_sapiens", gene_id=None, **retry):
    """``(gene id, every transcript id)`` of a gene (all biotypes), for use as background.

    By ``gene_id`` when given, and otherwise by the symbol ``gene``.  A symbol can name
    more than one gene -- ``lookup/symbol/CD99`` gives the chrY copy of this
    pseudoautosomal gene, whose sequence is the chrX copy's -- which is why a config
    records the gene it was proposed from.  Kept separate from
    :mod:`isoform_dominance.annotate`, which restricts itself to protein-coding
    transcripts: for identifiability the non-coding, retained-intron and NMD transcripts
    matter, because Salmon indexes them too.
    """
    seg = [urllib.parse.quote(str(x), safe="")      # a symbol can hold a space or a slash
           for x in ((gene_id.split(".")[0],) if gene_id else (species, gene))]
    path = ("/lookup/id/%s?expand=1" % tuple(seg) if gene_id
            else "/lookup/symbol/%s/%s?expand=1" % tuple(seg))
    info = ensembl.get_json(path, **retry)
    return info.get("id"), [t["id"].split(".")[0] for t in info.get("Transcript", [])]


def fetch_gene_transcript_ids(gene, species="homo_sapiens", gene_id=None, **retry):
    """Every transcript id of the gene :func:`fetch_gene` fetches."""
    return fetch_gene(gene, species, gene_id=gene_id, **retry)[1]


# --------------------------------------------------------------------------- #
# k-mer primitives
# --------------------------------------------------------------------------- #
def revcomp(seq):
    """Reverse complement of ``seq`` (upper-cased, N-safe)."""
    return seq.upper().translate(_COMPLEMENT)[::-1]


def canonical_kmer(km):
    """The lexicographically smaller of a k-mer and its reverse complement.

    Salmon's index is built on canonical k-mers, and an unstranded library gives
    fragments from either strand, so two sequences that differ only by orientation are
    *not* distinguishable in practice.  Canonicalising here keeps the uniqueness
    calculation from claiming a separability the quantifier does not have.
    """
    rc = revcomp(km)
    return km if km <= rc else rc


def kmers(seq, k, canonical=True):
    """Set of length-``k`` substrings of ``seq``.

    With ``canonical=True`` (the default) each k-mer is folded to the smaller of
    itself and its reverse complement.  Pass ``canonical=False`` for the strand-aware
    behaviour of releases up to v2.1.1.
    """
    return set(kmer_track(seq, k, canonical))


def kmer_track(seq, k, canonical=True):
    """Ordered list of the k-mers of ``seq``, one per start position.

    Unlike :func:`kmers` this preserves position, which is what the coverage,
    block-structure and read-model calculations need.
    """
    seq = seq.upper()
    n = len(seq)
    if k <= 0 or n < k:
        return []
    if canonical:
        # :func:`canonical_kmer` of each window, with one reverse complement of the whole
        # sequence: the window at i reads, on the other strand, as its slice at n - i - k
        rc = seq.translate(_COMPLEMENT)[::-1]
        return [min(seq[i:i + k], rc[n - i - k:n - i]) for i in range(n - k + 1)]
    return [seq[i:i + k] for i in range(n - k + 1)]


#: Length of the seeds the background-FASTA scan looks up before it reads a window.
FASTA_SEED = 16

#: A background-FASTA record longer than this, in nt, looks like a genome sequence rather
#: than a transcript: the longest human transcripts are around 0.1 Mb.
LONG_RECORD = 1_000_000


def _scan_records(path, query_kmers, k, canonical=True, exclude_ids=(), identical=None,
                  identical_out=None, seed=FASTA_SEED, decoys=(), stats=None):
    """Yield ``(record id, sequence, windows)`` for each record of a background FASTA.

    The sequence is the record's lines, stripped, joined and upper-cased; ``windows`` is
    the set of ``query_kmers`` the record holds, in the form ``query_kmers`` gives them.
    Records that :func:`scan_background_fasta` skips -- excluded, identical to a
    configured transcript, or with no sequence line -- are not yielded, and nothing is
    yielded for an empty query.  See :func:`scan_background_fasta` for the arguments and
    for how a window is found.

    ``decoys`` are record names, as Salmon's ``decoys.txt`` gives them: the header's
    first whitespace-delimited word, or that word's first ``|``-delimited field.  Such a
    record is skipped unread.  ``stats``, a dict, receives ``decoys_skipped`` (records),
    ``decoys_found`` (the set of names met) and ``long_records``, ``[(id, length)]`` of
    the records read that are longer than :data:`LONG_RECORD`.
    """
    decoys = frozenset(decoys)
    if stats is not None:
        stats.update(decoys_skipped=0, decoys_found=set(), long_records=[])
    query = set(query_kmers)
    if not query:
        return
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

    def _consume(chunks, tid):
        if not chunks or tid in exclude:
            return None
        seq = "".join(chunks).upper()
        if len(seq) in lengths and seq in identical:
            if identical_out is not None:
                identical_out[tid] = identical[seq]
            return None
        if stats is not None and len(seq) > LONG_RECORD:
            stats["long_records"].append((tid, len(seq)))
        return tid, seq, _hits(seq)

    chunks, tid, skip = [], None, False
    with io.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                rec = None if skip else _consume(chunks, tid)
                if rec is not None:
                    yield rec
                head = line[1:].strip()
                tid = head.split("|")[0].split()[0].split(".")[0] if head else None
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
        rec = None if skip else _consume(chunks, tid)
        if rec is not None:
            yield rec


def scan_background_fasta(path, query_kmers, k, canonical=True, exclude_ids=(),
                          identical=None, identical_out=None, *, seed=FASTA_SEED):
    """Return the subset of ``query_kmers`` that also occurs in a background FASTA.

    Streams the file and never materialises the background's own k-mer set, so a
    whole-transcriptome FASTA (GENCODE, or the transcript FASTA a Salmon index was built from)
    can be used as background; one record is held at a time, so memory is set by the
    query and the longest record, not the file (a single 20 Mb record: about 125 MB).
    Handles plain or gzipped input, told apart by the gzip magic bytes; ``exclude_ids`` drops records whose first
    ``|``- or whitespace-delimited field matches (version suffix ignored), which is
    how the transcripts under test are kept out of their own background.

    ``identical`` maps sequences, upper-cased and stripped, to the configured transcript
    each is: a record whose whole sequence is one of them is skipped too, and
    ``identical_out``, a dict, receives record id -> that transcript.  Salmon's index does
    the same with such a record: unless it is built with ``--keepDuplicates`` it keeps
    only the first of identical sequences, so the record competes with nothing.

    The answer is exact, and what folding every window of every record to its canonical
    form and looking it up gives; only the work differs.  Each query window is looked up
    in both orientations, so a record's windows are never folded.  And a record's window
    is read only where a ``seed``-long substring of it is one of the query's: with
    ``s = min(seed, k)`` and ``t = k - s + 1``, every window of length ``k`` holds the
    ``s``-mer that starts at the one multiple of ``t`` among its first ``t`` positions, so
    only the ``s``-mers at multiples of ``t`` are looked up, and the ``t`` windows that
    hold one that matches are then read whole.  ``seed`` changes the speed, never the
    answer.
    """
    seen = set()
    for _, _, hits in _scan_records(path, query_kmers, k, canonical, exclude_ids, identical,
                                    identical_out, seed):
        seen |= hits
    return seen


def scan_fasta_competitors(path, query_kmers, k, canonical=True, exclude_ids=(),
                           identical=None, identical_out=None, *, keep_ids=(),
                           keep_sequences=(), seed=FASTA_SEED, decoys=(), stats=None):
    """The records of a background FASTA that share a window with ``query_kmers``.

    ``{record id: (sequence, windows)}`` in file order, where ``windows`` is the set of
    query windows the record holds; :func:`scan_background_fasta`'s answer is the union
    of them.  A record that shares no window is left out unless its id (without its
    version) is in ``keep_ids`` or its sequence, upper-cased and stripped, is in
    ``keep_sequences``.  Only the records kept are held in memory.  A second
    record with an id already seen is left out when its sequence is the first's, and kept
    as ``<id>#2`` (``#3``, ...) when it is not; a record whose header gives no id is
    ``record<N>``, the N-th record scanned.  ``decoys`` and ``stats`` are described at
    :func:`_scan_records`; the other arguments are :func:`scan_background_fasta`'s.
    """
    keep = {str(i).split(".")[0] for i in keep_ids}
    keep_seqs = set(keep_sequences)
    keep_lengths = {len(x) for x in keep_seqs}
    out, n = {}, 0
    for tid, seq, hits in _scan_records(path, query_kmers, k, canonical, exclude_ids,
                                        identical, identical_out, seed, decoys, stats):
        n += 1
        if not (hits or tid in keep or (len(seq) in keep_lengths and seq in keep_seqs)):
            continue
        rid = tid if tid is not None else "record%d" % n
        if rid in out:
            if out[rid][0] == seq:
                continue
            m = 2
            while "%s#%d" % (rid, m) in out:
                m += 1
            rid = "%s#%d" % (rid, m)
        out[rid] = (seq, hits)
    return out


# --------------------------------------------------------------------------- #
# block structure
# --------------------------------------------------------------------------- #
def _blocks(flags):
    """Run-length encode a boolean array into (start, length) runs of True."""
    out, start = [], None
    for i, v in enumerate(flags):
        if v and start is None:
            start = i
        elif not v and start is not None:
            out.append((start, i - start))
            start = None
    if start is not None:
        out.append((start, len(flags) - start))
    return out


def coverage_stats(unique_flags, k):
    """Positional statistics for a transcript's group-unique k-mer starts.

    ``unique_flags[i]`` is True when the k-mer starting at position ``i`` is unique to
    the group.  Reported in *bases*, not k-mer counts: a k-mer count is hard to reason
    about (it depends on k and on how the unique region abuts shared sequence), while
    "how many bases of this transcript are uniquely attributable" is directly
    interpretable and is what the read model consumes.

    A base is uniquely attributable when *some* unique k-mer covers it, so the answer is
    the size of the **union** of the spans, not the sum of their lengths.  Two runs of
    unique starts separated by a gap of fewer than ``k`` positions have overlapping or
    touching spans: at ``k = 3`` with unique starts at 0 and 2 the spans are bases 0-2
    and 2-4, five bases, not six.  Summing run lengths double-counts the overlap and can
    drive ``unique_fraction`` above 1.  A ``block`` is likewise a maximal run of
    contiguous unique *bases*, which is what a read has to sit on, rather than a maximal
    run of unique k-mer starts.
    """
    flags = list(unique_flags)
    runs = _blocks(flags)
    # a run of r unique k-mer starts spans [s, s + r + k - 1); merge spans that touch
    # or overlap, because the bases under them are one contiguous unique stretch
    spans = []
    for s, r in runs:
        lo, hi = s, s + r + k - 1
        if spans and lo <= spans[-1][1]:
            spans[-1][1] = max(spans[-1][1], hi)
        else:
            spans.append([lo, hi])
    covered = sum(hi - lo for lo, hi in spans)
    block_lengths = sorted((hi - lo for lo, hi in spans), reverse=True)
    runs = spans
    return {
        "unique_length": int(covered),
        "n_blocks": len(runs),
        "max_block_length": int(block_lengths[0]) if block_lengths else 0,
        "block_lengths": [int(b) for b in block_lengths[:20]],
    }


# --------------------------------------------------------------------------- #
# read / fragment model
# --------------------------------------------------------------------------- #
def _fragment_length_grid(frag_mean, frag_sd, n_points=15, lo=None, hi=None):
    """Discretise a truncated normal fragment-length distribution to (length, weight)."""
    if frag_sd <= 0:
        return [(int(round(frag_mean)), 1.0)]
    lo = max(1.0, frag_mean - 3.0 * frag_sd) if lo is None else lo
    hi = frag_mean + 3.0 * frag_sd if hi is None else hi
    xs = np.linspace(lo, hi, n_points)
    w = np.exp(-0.5 * ((xs - frag_mean) / frag_sd) ** 2)
    w = w / w.sum()
    return [(int(round(x)), float(p))
            for x, p in zip(xs, w, strict=True) if p > 1e-6]


def informative_fraction(unique_flags, transcript_length, k,
                         read_length=DEFAULT_READ_LENGTH,
                         frag_mean=DEFAULT_FRAG_MEAN, frag_sd=DEFAULT_FRAG_SD,
                         paired=True):
    """P(a fragment from this transcript is unambiguously assignable to its group).

    A fragment is informative only if **a sequenced end** covers a whole group-unique
    k-mer.  Modelling the fragment rather than the reads would be wrong and optimistic:
    with a 200 nt fragment and 100 nt reads a unique stretch sitting in the middle of
    the fragment is never observed.  For a fragment ``[p, p+f)`` the observed windows
    are ``[p, p+L)`` and, when paired, ``[p+f-L, p+f)``.

    Evaluated exactly over all start positions by prefix sums, then averaged over a
    truncated-normal fragment-length distribution.
    """
    n_starts = len(unique_flags)
    if n_starts <= 0 or transcript_length < k:
        return 0.0
    flags = np.asarray(unique_flags, dtype=np.int64)
    csum = np.concatenate([[0], np.cumsum(flags)])

    def any_unique(a, b):
        """True where at least one unique k-mer starts in [a, b] (vectorised, inclusive)."""
        a = np.clip(a, 0, n_starts)
        b = np.clip(b + 1, 0, n_starts)
        return (csum[b] - csum[a]) > 0

    total = 0.0
    for f, w in _fragment_length_grid(frag_mean, frag_sd):
        f = min(f, transcript_length)
        if f < k:
            continue
        n_pos = transcript_length - f + 1
        if n_pos <= 0:
            continue
        p = np.arange(n_pos)
        L = min(read_length, f)
        hit = any_unique(p, p + L - k)
        if paired:
            back = p + f - L
            hit = hit | any_unique(back, back + L - k)
        total += w * float(hit.mean())
    return float(total)


def expected_informative_reads(informative_frac, tpm, transcript_length,
                               depth=DEFAULT_DEPTH, mean_efflen=DEFAULT_MEAN_EFFLEN,
                               frag_mean=DEFAULT_FRAG_MEAN):
    """Expected number of unambiguously assignable fragments for a class.

    ``TPM`` is normalised per molecule, so the fragment count scales with
    ``TPM x effective length``.  With ``sum(TPM) = 1e6`` over the transcriptome, the
    expected fragments from a transcript are approximately
    ``depth x TPM x efflen / (1e6 x mean_efflen)``; ``mean_efflen`` is the
    TPM-weighted mean effective length of the library and is exposed because it is the
    one term that cannot be derived from the gene alone.
    """
    efflen = effective_length(transcript_length, frag_mean)
    reads = depth * (tpm * efflen) / (1e6 * max(1.0, mean_efflen))
    return float(reads * informative_frac)


def effective_length(transcript_length, frag_mean=DEFAULT_FRAG_MEAN):
    """``max(1, L - frag_mean + 1)``: the start positions a fragment of mean length has."""
    return max(1.0, transcript_length - frag_mean + 1.0)


#: Smallest ``|log2_efflen_ratio|`` for which the CLI states a skew direction.  In the
#: 49-gene simulation the direction was scored on the genes whose class effective lengths
#: differ by more than 1.23x (log2 0.2987); 0.3 selects the same 39 genes there.  Below it
#: nothing was measured, so nothing is said.
DIRECTION_MIN_ABS_LOG2_EFFLEN_RATIO = 0.3


def class_efflen_ratio(lengths_a, lengths_b, frag_mean=DEFAULT_FRAG_MEAN):
    """``(mean_efflen_a, mean_efflen_b, log2(mean_efflen_a / mean_efflen_b))``.

    The mean is the plain mean over each class's transcripts of
    :func:`effective_length`.  ``None`` for all three when either class is empty.

    In a 49-gene simulation (Salmon, monotone positional coverage skew) the quantifier
    split the ambiguous mass between the classes by effective length, so this ratio
    tracked the size of the error and, with the direction of the skew, its sign: a
    5'-skewed library tended to inflate the shorter class, a 3'-skewed one the longer.
    No other quantifier, panel or form of skew was tested.  The tendency assumes the
    classes are ambiguous where the skew concentrates reads; a shorter class distinguished
    at that end is an untested interaction under which it can reverse -- LEPR, shorter
    class distinguished at its 5' end, ran against it under 5 of 6 skews.  The mechanism
    is about each transcript's effective length, so the transcript count must not enter.
    Five summaries were scored against realised bias on that simulation (Spearman at the
    two extreme skews, and how often the 5' sign came out right among genes whose ratio
    exceeds 1.23x):

    ===================  ==========  ==========  =====
    summary              b = +2.0    b = -2.0    sign
    ===================  ==========  ==========  =====
    plain mean           -0.547      +0.169      35/39
    harmonic mean        -0.541      +0.164      36/41
    class total          -0.369      +0.285      35/47
    minimum efflen       -0.398      +0.101      32/39
    sum of 1/efflen      -0.369      -0.058      24/45
    ===================  ==========  ==========  =====

    The plain mean is best or tied.  The class total, which weights by transcript count,
    is worse at the extremes and reads -0.344 at b = 0, where the answer has to be null:
    it is picking up class size.
    """
    if not lengths_a or not lengths_b:
        return None, None, None
    mean_a = sum(effective_length(n, frag_mean) for n in lengths_a) / len(lengths_a)
    mean_b = sum(effective_length(n, frag_mean) for n in lengths_b) / len(lengths_b)
    return float(mean_a), float(mean_b), math.log2(mean_a / mean_b)


def window_positions(flags):
    """Where the flagged windows start, as a fraction of the transcript's own length.

    ``flags[i]`` is True when the window starting at ``i`` distinguishes the class.  The
    result runs from 0 (the window at the 5' end) to 1 (the last window, at the 3' end):
    ``i / (len(flags) - 1)``.

    Each transcript is measured in its own coordinates, not the gene's.  Positional
    coverage skew acts on each transcript separately -- reads pile up toward each
    molecule's own end -- so that is the coordinate system the mechanism runs in.  A
    window 500 nt from the 3' end of a 700-nt transcript and one 500 nt from the 3' end of
    a 3.2-kb transcript are in very different places as far as the skew is concerned.
    """
    last = max(1, len(flags) - 1)
    return [i / last for i, f in enumerate(flags) if f]


def position_summary(positions):
    """``{n, q1, median, q3}`` of :func:`window_positions` pooled over a class.

    ``n`` counts positions, not distinct windows: a distinguishing window carried by two
    of the class's transcripts has a position in each, and counts twice.
    """
    if not positions:
        return {"n": 0, "q1": None, "median": None, "q3": None}
    q1, med, q3 = np.percentile(positions, [25, 50, 75])
    return {"n": len(positions), "q1": float(q1), "median": float(med), "q3": float(q3)}


def counting_noise_floor(n_informative_a, n_informative_b, n_donors=1):
    """Counting-noise standard error of the per-sample log2 class ratio.

    Poisson counting error alone; biological and technical variation add on top, so
    this is a floor on the achievable precision, never an estimate of it.  Returned
    as NaN when either class has no informative reads, because the ratio is then
    undefined rather than merely imprecise.
    """
    if n_informative_a <= 0 or n_informative_b <= 0:
        return {"log2_ratio_se": float("nan"), "min_resolvable_log2fc": float("nan")}
    se = math.sqrt(1.0 / n_informative_a + 1.0 / n_informative_b) / math.log(2.0)
    n_donors = max(1, int(n_donors))
    return {
        "log2_ratio_se": se,
        "min_resolvable_log2fc": 1.96 * se / math.sqrt(n_donors),
    }


def class_coherence(tracks, ids):
    """How much the transcripts assigned to one class look alike, as window Jaccard.

    ``annotate`` groups protein-coding transcripts by their 3' terminal-exon splice
    acceptor and its docstring asks the user to review the proposal.  This is the number
    to review it with: a class whose members share almost no windows is not a functional
    group that happens to be hard to measure, it is two unrelated transcripts that share
    one acceptor, and every precision figure computed for it describes a quantity nobody
    wants.  Across a 49-gene survey three classes came in under 0.05 (TPI1 0.004, whose
    two members are 374 nt and 2217 nt long; CD44 0.009; DMD 0.038) while the median
    class sat near 0.7.

    Low coherence is a reason to revisit the grouping, not a verdict: a class can be
    perfectly coherent and still be hard to measure, which is what most of that survey's
    difficult genes turned out to be.
    """
    members = [t for t in ids if t in tracks]
    if len(members) < 2:
        return {"median_jaccard": None, "min_jaccard": None, "n_pairs": 0}
    sets = {t: frozenset(tracks[t]) for t in members}
    js = []
    for i, a in enumerate(members):
        for b in members[i + 1:]:
            union = len(sets[a] | sets[b])
            js.append(len(sets[a] & sets[b]) / union if union else 1.0)
    js.sort()
    mid = len(js) // 2
    median = js[mid] if len(js) % 2 else 0.5 * (js[mid - 1] + js[mid])
    return {"median_jaccard": float(median), "min_jaccard": float(js[0]),
            "n_pairs": len(js)}


def gls_covariance(A, theta):
    """Covariance of the GLS estimate of ``theta`` under Poisson counts: ``(A'WA)^+``.

    :func:`estimability` reports ``sqrt(c'(A'A)^+ c)``, which is the standard-deviation
    multiplier of the BLUE only under ``Var(y) = sigma^2 I``.  Fragment counts are not
    homoskedastic: ``Var(y_c) = E[y_c]``, so the weight is ``W = diag(1/mean)``.

    The difference is not cosmetic.  Over 41 estimable genes the homoskedastic factor had
    a median of 65.0 and a range of 850x; the same genes under Poisson weighting had a
    median of 16.4 and a range of 157x.  The *ranking* barely moved (Spearman 0.93), so
    this does not rescue a fixed threshold on the factor -- what it does is put the number
    on a scale that can be turned into a detectable effect.

    ``theta`` is the vector of expected fragment counts per transcript.  It is what the
    analysis is trying to estimate, so a plug-in is unavoidable; :func:`analyze` passes a
    flat one from the stated design.  Real genes are skewed toward one isoform and the
    weights move with that, in a direction this function does not attempt to predict.
    """
    A = np.asarray(A, dtype=float)
    theta = np.asarray(theta, dtype=float)
    if A.size == 0:
        return None
    mean = A @ theta
    keep = mean > 1e-12
    if not keep.any():
        return None
    aw = A[keep] / np.sqrt(mean[keep])[:, None]
    return np.linalg.pinv(aw.T @ aw, rcond=1e-10)


def gls_relative_se(cov, c, theta):
    """Standard error of ``c'theta``, relative to it -- the SE of ``log(c'theta)``."""
    if cov is None or not np.any(c):
        return float("inf")
    var = float(c @ cov @ c)
    est = float(c @ theta)
    if var < 0.0 or not est:
        return float("inf")
    return math.sqrt(var) / abs(est)


def log_ratio_se(cov, c_a, c_b, theta):
    """Standard error of ``log(A/B)`` for two class totals, by the delta method.

    The contrast's estimand is the *ratio* of the class totals, so its precision is not
    ``SD(A-B)`` rescaled: a difference can sit near zero for reasons that have nothing to
    do with how well either total is determined, and the two totals are correlated because
    they are solved from one system.  With ``Sigma`` the GLS covariance,

        Var(log A/B) = c_a'Sigma c_a / A^2 + c_b'Sigma c_b / B^2 - 2 c_a'Sigma c_b / (A B)

    which is what makes the returned figure comparable to the per-class ones and
    convertible into a fold change.
    """
    if cov is None:
        return float("inf")
    a = float(c_a @ theta)
    b = float(c_b @ theta)
    if not a or not b:
        return float("inf")
    va = float(c_a @ cov @ c_a) / (a * a)
    vb = float(c_b @ cov @ c_b) / (b * b)
    cab = float(c_a @ cov @ c_b) / (a * b)
    var = va + vb - 2.0 * cab
    if var < 0.0:                       # numerically negative only when it is ~0
        var = 0.0
    return math.sqrt(var)


LINEARISATION_LIMIT = 0.3
"""Above this standard error on the log scale the first-order figures stop being values.

Checked against simulation -- Poisson counts, GLS fit, the sample SD of ``log2(A/B)`` over
4000 draws.  Below the limit the delta-method figure is the answer: PIK3CA came back at
1.002 and 0.992 times the simulated SD at 1e4 and 1e6 fragments, FLT1 at 1.034.  Above it
the two part company -- NR1H3 1.244, TPI1 0.540 -- because the expansion is first order and
because the unconstrained GLS total can go negative there, which truncates the simulation
too.  Neither number is trustworthy in that regime; both still say the comparison is out of
reach, which is the only thing being claimed.
"""


def min_resolvable_log2fc(relative_se, n_donors=1):
    """Smallest |log2 fold change| a 95% interval excludes zero for, at this design.

    The companion to :func:`counting_noise_floor`, which answers the same question for an
    estimator that counts only unambiguously assignable reads.  This one is for the
    estimator that inverts the whole system, which is what ``conditioning_factor``
    describes -- the two numbers in a report have always belonged to different estimators.

    Still a floor: Poisson counting error only, biological and technical variation on top.
    And a bound on spread, not on accuracy: it assumes uniform coverage and a correctly
    specified compatibility model.  Under positional coverage skew the estimate can be
    biased well past it while the replicate SD stays below the predicted SE, so replicate
    agreement does not reveal the bias; the README gives the simulation numbers.
    A figure past :data:`LINEARISATION_LIMIT` on the log scale reads as "not resolvable at
    this design", not as a calibrated value; ``analyze`` flags those as ``beyond_linear``.

    Not comparable to :func:`counting_noise_floor`'s figure, which is the precision of
    ``log2(n_a/n_b)`` for the informative-read counts of each class's *best single
    transcript*.  That ratio equals the class ratio only when the two classes share an
    informative fraction, so the older number can look much better than the class
    comparison is -- on TPI1 it reads 0.76 where the class totals themselves are not
    resolvable at all.
    """
    if not math.isfinite(relative_se):
        return float("inf")
    return 1.96 * relative_se / (math.log(2.0) * math.sqrt(max(1, int(n_donors))))


# --------------------------------------------------------------------------- #
# estimability of the class-collapsed system
# --------------------------------------------------------------------------- #
def compatibility_matrix(tracks, transcript_ids):
    """Design matrix of the fragment-class system, one row per compatibility class.

    ``tracks`` maps transcript id to its ordered window list (k-mers, or longer
    windows -- see ``window`` in :func:`analyze`).  Windows sharing a *signature* --
    the set of transcripts that contain them -- are indistinguishable to a quantifier
    and are collapsed into one row.  ``A[c, t]`` is the probability that a window
    drawn uniformly from transcript ``t`` falls in class ``c``, so
    ``E[count_c] = sum_t theta_t A[c, t]`` up to a shared depth factor.  Each column
    therefore sums to 1 by construction.

    That probability is over *positions*, not over distinct window sequences.  A window
    that occurs twice in a transcript is drawn twice as often, and repeated windows are
    the rule rather than the exception in cDNA -- A-rich 3' ends, tandem repeats, Alu
    elements in long UTRs.  Counting distinct sequences instead under-weights exactly
    the shared classes that absorb the most fragments.

    **This is a surrogate, not the observation model of a sequencing run.**  What an
    actual paired-end library observes is decided by read length, the fragment-length
    distribution, and the fact that only the two ends of a fragment are sequenced.  No
    claim is made relating the two systems -- not that a verdict here transfers to
    reads, and not that it is conservative in either direction.  Several attempts to
    state such a relation failed, and the counterexample that ended the last of them is
    ``test_a_longer_window_is_a_different_system_not_a_sharper_one``: four transcripts
    whose class contrast is estimable at window 3 and 4, not at 5 and 6, and estimable
    again at 7 and 8.  Raising ``window`` gives a different system, not a sharper one.
    Read a verdict from this matrix as a statement about a sequence-derived screening
    surrogate at a stated ``window``, never as a property of the data; whether it
    predicts what a quantifier recovers is the empirical question the simulation study
    addresses.

    One caveat on the column sums: a transcript shorter than ``window`` has no windows
    and gets an all-zero column, so ``A`` is column-stochastic only when every
    transcript is at least ``window`` long.
    """
    # A window's signature, built without a container for the windows only one transcript
    # has -- most of them, and millions once many background records join -- and with one
    # frozenset per signature, shared by all its windows, which a dictionary lookup then
    # tells apart by identity rather than element by element.
    sig, shared = {}, []
    for tid in transcript_ids:
        for w in set(tracks.get(tid, [])):
            have = sig.get(w)
            if have is None:
                sig[w] = tid                    # this transcript's alone, so far
            elif isinstance(have, list):
                have.append(tid)
            else:
                sig[w] = [have, tid]
                shared.append(w)
    interned = {}
    for w in shared:
        key = frozenset(sig[w])
        sig[w] = interned.setdefault(key, key)

    # positions carrying each signature, per transcript -- see the docstring on why
    # this is not a count of distinct window sequences
    pos_counts = {}
    for tid in transcript_ids:
        for key, n in Counter(map(sig.__getitem__, tracks.get(tid, []))).items():
            if not isinstance(key, frozenset):  # a window no other transcript has
                key = frozenset((key,))
                key = interned.setdefault(key, key)
            pos_counts[(key, tid)] = pos_counts.get((key, tid), 0) + n
    class_keys = set(interned)

    idx = {t: j for j, t in enumerate(transcript_ids)}
    n_windows = {t: max(1, len(tracks.get(t, []))) for t in transcript_ids}
    classes = sorted(class_keys, key=lambda s: (-len(s), sorted(s)))
    A = np.zeros((len(classes), len(transcript_ids)), dtype=float)
    for i, key in enumerate(classes):
        for t in key:
            A[i, idx[t]] = pos_counts.get((key, t), 0) / n_windows[t]
    return A, [sorted(c) for c in classes]


def estimability(A, c, rcond=1e-10, *, svd=None):
    """Is the linear functional ``c'theta`` estimable from ``E[y] = A theta``, and how well?

    ``svd`` is ``np.linalg.svd(A.T, full_matrices=False)[:2]``, for a caller that asks
    about several functionals of one ``A``: the answer is the one this function gives
    without it, computed once.

    **Structure.**  ``c'theta`` is estimable exactly when ``c`` lies in the row space of
    ``A``.  The residual of projecting ``c`` onto that row space, relative to ``||c||``,
    is reported as ``residual``.  This part is exact and model-free: it depends only on
    which fragment classes are observable, not on how their counts are distributed.

    **Conditioning.**  ``conditioning_factor`` is ``sqrt(c' (A'A)^+ c)``.  Read it as a
    *structural conditioning proxy*, not as the standard error of anything.  Note the
    square root: the variance factor is ``c'(A'A)^+ c`` and this is its
    standard-deviation counterpart, so that ``SD(BLUE of c'theta) = sigma *
    conditioning_factor`` would hold under ``Var(y) = sigma^2 I`` -- an assumption this
    package does not make and that a short-read quantifier does not satisfy: fragment
    counts are heteroskedastic, and
    Salmon's rich equivalence classes carry per-transcript conditional probabilities and
    bias weights rather than the 0/1-derived compatibilities used to build ``A`` here.
    What the quantity does capture is the geometry: a contrast direction that is nearly
    degenerate in the observable classes has a large factor and will be poorly
    determined however the noise is actually distributed, which is precisely the case a
    rank test alone waves through.

    Whether the proxy predicts realised quantification error is an empirical question,
    not a theorem, and is the subject of the simulation study tracked in the repository
    issues.  Until that is answered, treat it as a screening quantity and do not report
    it as a variance.
    """
    A = np.asarray(A, dtype=float)
    c = np.asarray(c, dtype=float)
    if A.size == 0:
        # no observable classes: the row space is {0}, so only the zero functional is
        # estimable, and it is estimated by the constant 0 with no error
        zero = not np.any(c)
        return {"estimable": zero, "residual": 0.0 if zero else 1.0,
                "conditioning_factor": 0.0 if zero else float("inf"), "rank": 0}

    # Row space of A == column space of A.T, so the columns of ``u`` are the right
    # singular vectors of A and ``s`` its singular values.  Rank and conditioning are
    # computed from the SAME truncation: taking the rank from sigma(A) while taking the
    # conditioning from pinv(A'A) -- whose rcond is relative to sigma(A)^2 -- puts a
    # band of directions on both sides of the line at once, where a contrast is
    # declared estimable and its conditioning direction is simultaneously projected
    # away, reporting 0.0 (the best possible score) for the worst-conditioned case.
    u, s = svd if svd is not None else np.linalg.svd(A.T, full_matrices=False)[:2]
    tol = max(A.shape) * (s[0] if s.size else 0.0) * rcond
    r = int((s > tol).sum())
    if not np.any(c):
        return {"estimable": True, "residual": 0.0,
                "conditioning_factor": 0.0, "rank": r}
    basis = u[:, :r]
    coef = basis.T @ c                      # coordinates of c in the retained row space
    proj = basis @ coef
    denom = np.linalg.norm(c) or 1.0
    residual = float(np.linalg.norm(c - proj) / denom)
    estimable = residual < 1e-8
    if estimable and r:
        # c'(A'A)^+ c = sum_i (v_i'c)^2 / sigma_i^2 over the retained directions
        cond = float(math.sqrt(float(np.sum((coef[:r] / s[:r]) ** 2))))
    else:
        # a contrast with a component outside the row space has no unbiased estimator,
        # so there is no finite factor to report for it
        cond = float("inf")
    return {
        "estimable": estimable,
        "residual": residual,
        "conditioning_factor": cond,
        "rank": r,
    }


# --------------------------------------------------------------------------- #
# top-level analysis
# --------------------------------------------------------------------------- #
def _verdict(entry, tau, min_reads=None, min_log2fc=None):
    """Grade one estimand, and say why.

    ``entry`` may be a class total or the contrast; ``entry["label"]`` names which, so
    the reason string does not assert the wrong one.  Structure first: a functional
    outside the row space of *this surrogate system*
    cannot be recovered from it at any depth -- which is a statement about the
    surrogate, not about the data; see :func:`compatibility_matrix`.  Then precision,
    from two independent directions -- an ill-conditioned contrast (a large conditioning
    factor, which is a standard-deviation factor and not a variance) and simple lack of
    informative fragments.  The
    second is why "has at least one unique k-mer" is not a usable gate: a class whose
    uniqueness is thirty junction k-mers is structurally estimable and practically
    hopeless.
    """
    reasons = []
    if not entry.get("estimable"):
        reasons.append("%s outside the row space of the compatibility system"
                       % entry.get("label", "estimand"))
        return "not_identifiable", reasons
    if min_log2fc is not None:
        # An effect size the user can defend, in place of a conditioning number nobody
        # can: tau has no calibrated value -- over a 49-gene survey the median gene sat
        # at 65.0 against a default tau of 10.0, so the default rejects 85% of what it
        # is applied to and any other fixed value simply sorts genes by how many
        # transcripts they have annotated.  See docs and the CHANGELOG entry.
        got = entry.get("min_resolvable_log2fc")
        if got is None or not math.isfinite(got):
            reasons.append("%s has no finite resolvable effect size at this design"
                           % entry.get("label", "estimand"))
        elif got > min_log2fc:
            reasons.append("smallest resolvable |log2FC| %.2f exceeds the requested %.2f"
                           % (got, min_log2fc))
        elif entry.get("beyond_linear"):
            # the first-order figure is not a value past the limit, so it cannot resolve
            # anything however small it reads
            reasons.append("%s is beyond the linearisation limit (relative SE %.2f > %.1f), "
                           "so not resolvable at this design"
                           % (entry.get("label", "estimand"), entry["gls_relative_se"],
                              LINEARISATION_LIMIT))
    elif entry.get("conditioning_factor", 0.0) > tau:
        reasons.append("conditioning factor %.1f exceeds tau=%.1f"
                       % (entry["conditioning_factor"], tau))
    n = entry.get("expected_informative_reads")
    if min_reads is not None and n is not None and n < min_reads:
        reasons.append("expected informative reads %.0f below %.0f" % (n, min_reads))
    return ("weakly_identifiable" if reasons else "identifiable"), reasons


def analyze(config, k=DEFAULT_K, sequences=None, *,
            canonical=True, window=None,
            background_sequences=None, background_fasta=None, decoys=None,
            background_gene_transcripts="auto", species=None,
            read_length=DEFAULT_READ_LENGTH, frag_mean=DEFAULT_FRAG_MEAN,
            frag_sd=DEFAULT_FRAG_SD, paired=True, depth=DEFAULT_DEPTH,
            mean_efflen=DEFAULT_MEAN_EFFLEN, tpm=DEFAULT_TPM, n_donors=1,
            conditioning_tau=DEFAULT_CONDITIONING_TAU,
            min_informative_reads=DEFAULT_MIN_INFORMATIVE_READS,
            min_log2fc=None, ensembl_release=None, inputs_out=None,
            keep_duplicates=False,
            retries=DEFAULT_RETRIES, retry_wait=DEFAULT_RETRY_WAIT):
    """Assess whether the configured isoform classes are measurable by short reads.

    Parameters
    ----------
    config
        The pipeline config dict; ``groups`` and ``primary_comparison`` are used, and
        ``ensembl_release``, when present, is carried into the report.
    k, window
        k-mer length, and the window length used to build the compatibility system
        (defaults to ``k``; a different window gives a different system, not a
        uniformly sharper one -- the rank is not monotone in it, see
        :func:`compatibility_matrix`)
    sequences
        ``{transcript_id: cdna}``.  Anything missing is fetched from Ensembl.
    background_sequences, background_fasta, background_gene_transcripts
        What uniqueness is judged against, beyond the other configured groups, and the
        columns of the compatibility system beside the configured transcripts.  Pass a
        FASTA -- the transcript FASTA the Salmon index was built from, without the genome
        decoys (or pass ``decoys``) -- for the honest
        whole-index answer.  Each of its records that shares a window with a configured
        transcript is a column, as a transcript of the gene background is, so a competitor
        gives one answer whether it comes from the FASTA or from ``background_sequences``.
        A record that shares no window with a configured transcript is not a column, and
        is not followed through the background, although it can bear on the configured
        estimands through a background transcript whose windows it shares: in one
        constructed case a configured transcript lying wholly inside a background
        transcript was estimable, and was not once a record held every window of the
        background transcript outside it.  A record with the id (without version) of a
        gene-background transcript is that transcript, once, with the FASTA's sequence
        when the two differ.  Salmon's default index keeps one of identical sequences, so
        a background sequence -- of the gene, from ``background_sequences`` or a FASTA
        record -- that is a configured transcript's is not counted, and one that is
        another background sequence's is counted once, as the first of them (the gene
        background in id order, then the FASTA in file order).  ``keep_duplicates=True``
        counts every copy, for an index built with ``--keepDuplicates``.
        ``background_gene_transcripts`` defaults to ``"auto"``:
        the gene's remaining transcripts are fetched and used when the caller is
        already relying on Ensembl for sequence, and skipped when sequences were
        supplied offline (so an offline call never blocks on the network).  ``True``
        forces the fetch, ``False`` restores the pre-v2.2 behaviour of comparing the
        configured groups only.  A gene symbol Ensembl does not know (HTTP 400/404)
        leaves the gene background empty; a ``gene_id`` it does not know is a
        ``ValueError``, and any other failure to fetch it is raised,
        because a background that is only partly fetched silently gives a different
        answer.
    decoys
        Salmon's ``decoys.txt`` -- one record name per line -- when ``background_fasta``
        is the gentrome of a decoy-aware index, transcripts and genome.  Those records are
        skipped unread, and their number reported.  Salmon sets aside only the reads that
        map better to a decoy than to any transcript, so a decoy competes with no read the
        transcripts explain as well: an exon's reads match the genome as well as the
        transcript, and judged against the genome every window inside an exon would lose
        its uniqueness.  Without ``decoys``, the records longer than :data:`LONG_RECORD`
        are reported as ``fasta_long_records``; pass the transcript FASTA the index was
        built from, or its ``decoys.txt``.
    read_length, frag_mean, frag_sd, paired, depth, mean_efflen, tpm, n_donors
        The sequencing design the report should be conditioned on.
    conditioning_tau, min_informative_reads
        Thresholds separating ``identifiable`` from ``weakly_identifiable``.
        ``conditioning_tau`` has no calibrated value and is kept as the default only for
        compatibility: a 49-gene survey put the median gene at a conditioning factor of
        65.0 against the default 10.0, and no fixed threshold does better than sorting
        genes by their annotated transcript count.  Prefer ``min_log2fc``.
    min_log2fc
        The smallest class-ratio change, as ``|log2 fold change|``, the caller needs to
        resolve.  When given it replaces ``conditioning_tau`` in the verdict, and the
        comparison is against :func:`min_resolvable_log2fc` computed from the GLS
        standard error of the *whole* system under Poisson counts -- not from
        :func:`counting_noise_floor`, which describes a unique-read-counting estimator
        and is reported alongside for what it is.
    ensembl_release
        The Ensembl release to fetch sequence and the gene background from, when
        anything has to be fetched.  None means the release ``rest.ensembl.org``
        currently serves; an earlier one is read from Ensembl's REST archive (see
        :func:`isoform_dominance.ensembl.resolve_server`).  Nothing is resolved, and no
        request made, when every sequence was supplied.
    inputs_out
        A dict to fill with the sequence this run used: ``sequences`` (the configured
        transcripts), ``background_sequences`` (the gene background, after the configured
        transcripts are removed from it), ``gene_background`` (whether this run had one, fetched
        or supplied), ``fetched_release``, ``sequence_sources`` (each id's ``"supplied"`` or
        ``"fetched:<release>"``), and the ``k``, ``window`` and ``canonical`` the system was
        built at -- none of which is in the config.
        Written to a file by :func:`isoform_dominance.io.save_inputs`, it repeats the run
        with no request at all, after the release's REST archive is gone.
    retries, retry_wait
        Each Ensembl request is retried up to ``retries`` times after the first
        attempt, waiting ``retry_wait`` seconds and doubling each time (an HTTP 429
        waits for its ``Retry-After``); defaults 5 and 1.0.  cDNA is fetched 50
        transcripts per request.  See :mod:`isoform_dominance.ensembl`.

    Returns
    -------
    dict
        ``groups`` (per-class sequence, read-model and estimability numbers),
        ``primary_comparison``, ``contrast`` (estimability of ``s_A - s_B``, with
        ``log2_efflen_ratio`` and ``class_mean_efflen`` from :func:`class_efflen_ratio`,
        ``efflen_direction_in_band``, and ``distinguishing_window_position`` per class
        from :func:`position_summary`),
        ``verdict`` with its ``reasons``, ``background`` (``gene_transcripts`` and
        ``n_background_transcripts``, the gene background's columns;
        ``fasta_competitors``, the FASTA records that are columns, each with the number of
        distinct windows it shares with the configured transcripts, and
        ``n_fasta_competitors``; ``fasta_sha256`` and ``decoys_sha256``, the SHA-256 of
        the files' bytes; ``sequence_from_fasta``, the gene-background transcripts
        whose column is the FASTA's sequence; ``identical_to_configured`` and
        ``identical_to_background``, the background sequences left out as copies, each
        mapped to the transcript it equals, and ``identical_source``, where each came
        from: ``"gene"``, ``"sequences"`` or ``"fasta"`` -- all three None with
        ``keep_duplicates``), ``annotation`` (the config's
        ``ensembl_release``, and ``fetched_release``, the release any sequence was
        fetched from in this run -- None for both when absent), ``gene_total``
        (estimability of the sum of every column, with the transcripts that have no
        window at all),
        ``effect_resolvable`` (whether both class totals and the contrast resolve
        ``min_log2fc``, none of them ``beyond_linear``; None without it), and -- for callers written against v2.1 --
        ``primary_distinguishable``.  Note that ``verdict`` supersedes
        ``primary_distinguishable``: a class with no unique k-mer of its own is still
        estimable when a class it is nested inside has unique sequence, and a class
        whose only unique sequence is a handful of junction k-mers is structurally
        estimable but practically not.
    """
    groups = config["groups"]
    if not groups:
        raise ValueError("config['groups'] is empty; nothing to test.")
    window = int(window or k)
    if decoys and not background_fasta:
        raise ValueError("--decoys names records of a --background-fasta, and none was given")

    pc = config.get("primary_comparison", list(groups)[:2])
    if len(pc) < 2:
        raise ValueError(
            "primary_comparison must name two isoform groups; got %r" % (pc,))
    missing = [g for g in pc if g not in groups]
    if missing:
        raise ValueError(
            "primary_comparison names group(s) not in config['groups']: %r" % (missing,))

    shared = io.shared_transcripts(groups)
    if shared:
        # in the two compared classes its +1 and -1 cancel and the contrast silently drops
        # it; in any two, its column enters the system twice
        t, gs = sorted(shared.items())[0]
        raise ValueError("transcript %s is in groups %s%s; a transcript belongs to one class"
                         % (t, " and ".join('"%s"' % g for g in gs),
                            " (and %d more transcript(s) are in two groups)"
                            % (len(shared) - 1) if len(shared) > 1 else ""))
    if window > read_length:
        # no read holds a whole window, so no fragment is informative and every class
        # reads as unmeasurable for a reason that is not the gene's
        raise ValueError("window %d exceeds the read length %d: no read can hold a whole "
                         "window, so no fragment would count as informative; pass a "
                         "--window no longer than --read-length" % (window, read_length))
    group_ids = {g: [t.split(".")[0] for t in ids] for g, ids in groups.items()}
    needed = [t for ids in group_ids.values() for t in ids]
    seqs = {tid.split(".")[0]: s for tid, s in (sequences or {}).items()}

    # ---- background transcripts of the same gene -------------------------- #
    bg_seqs = {t.split(".")[0]: s for t, s in (background_sequences or {}).items()}
    supplied = set(seqs) | set(bg_seqs)         # for inputs_out's sequence_sources
    if background_gene_transcripts == "auto":
        # only reach for the network when we are already going there for sequence
        background_gene_transcripts = any(t not in seqs for t in needed)
    net = {"retries": retries, "retry_wait": retry_wait}
    fetch_gene_background = bool(background_gene_transcripts and config.get("gene")
                                 and not bg_seqs)
    # The release sequence comes from need not be the one the config was annotated
    # against.  Record it whenever this run takes sequence from a server, and never
    # otherwise: a run on supplied sequence did not use one.
    if fetch_gene_background or any(t not in seqs for t in needed):
        net["server"] = ensembl.resolve_server(ensembl_release, **net)
        fetched_release = ensembl.fetch_release(**net)
    else:
        fetched_release = None
    # the configured transcripts first: when the release lacks them, that is what to say,
    # not that the gene background around them belongs to another gene
    missing = [t for t in needed if t not in seqs]
    if missing:
        got = ensembl.fetch_cdna_batch(missing, **net)
        absent = [t for t in missing if t not in got]
        if absent:
            raise ValueError("Ensembl release %s has no cDNA for %s, which the config names"
                             % (fetched_release, ", ".join(absent)))
        seqs.update(got)
    bg_gene_id = None
    if fetch_gene_background:
        try:
            bg_gene_id, all_ids = fetch_gene(
                config["gene"], species or config.get("species", "homo_sapiens"),
                gene_id=config.get("gene_id"), **net)
        except HTTPError as e:
            if e.code not in (400, 404):      # 400 is Ensembl's "no such symbol" or "id"
                raise
            if config.get("gene_id"):
                # the config names this gene outright; going on without its background
                # would judge the classes against nothing and call that a verdict
                raise ValueError("Ensembl release %s has no gene %s, which the config's "
                                 "gene_id names" % (fetched_release,
                                                    config["gene_id"].split(".")[0])) from e
            all_ids = []                      # a symbol Ensembl does not know
        if all_ids and not set(all_ids) & set(needed):
            # a symbol that names more than one gene gave another one: its transcripts are
            # not this gene's, and for a pseudoautosomal gene they are this gene's sequence
            raise ValueError(
                "Ensembl gives %s for %s, and none of its %d transcripts is one the config "
                "names, so the gene background would be another gene's. %s"
                % (bg_gene_id, "gene_id %s" % config["gene_id"] if config.get("gene_id")
                   else "the symbol %s" % config["gene"], len(all_ids),
                   "Check the config's gene_id." if config.get("gene_id") else
                   "A symbol can name more than one gene: add the config's \"gene_id\" "
                   "(`annotate` records it), or pass --no-gene-background."))
        # all or nothing: this used to swallow any error part-way through and carry
        # on with whatever had arrived, which is a different answer, silently
        bg_seqs.update(ensembl.fetch_cdna_batch(
            [t for t in all_ids if t not in needed], **net))
    bg_seqs = {t: s for t, s in bg_seqs.items() if t not in needed}
    if inputs_out is not None:
        inputs_out.update(sequences={t: seqs[t] for t in needed},
                          background_sequences=dict(bg_seqs),
                          # whether the rest of the gene was this run's background: a rerun
                          # under another grouping must not move transcripts into a
                          # background the run never had
                          gene_background=bool(fetch_gene_background
                                               or background_sequences is not None),
                          gene_id=bg_gene_id,
                          keep_duplicates=keep_duplicates,
                          fetched_release=fetched_release,
                          sequence_sources={
                              t: "supplied" if t in supplied else "fetched:%s" % fetched_release
                              for t in list(needed) + sorted(bg_seqs)},
                          k=k, window=window, canonical=canonical)

    # ---- window tracks ---------------------------------------------------- #
    tracks = {t: kmer_track(seqs[t], window, canonical) for t in needed}
    group_windows = {g: set().union(*(set(tracks[t]) for t in ids)) if ids else set()
                     for g, ids in group_ids.items()}

    # ---- an external FASTA background, streamed --------------------------- #
    # A record that shares a window with a configured transcript is a column of the system,
    # as the gene's other transcripts are.  Through 2.4.1 it took windows from the
    # uniqueness layer only, so the same competitor gave a smaller rank and a smaller
    # min |log2FC| passed as a FASTA than passed with background_sequences.
    source = "gene" if fetch_gene_background else "sequences"
    norm = {t: bg_seqs[t].strip().upper() for t in bg_seqs}
    # id -> (sequence, where it came from), the gene background first, sorted
    candidates = {t: (norm[t], source) for t in sorted(bg_seqs)}
    configured = {}
    for t in sorted(set(needed)):
        configured.setdefault(seqs[t].strip().upper(), t)
    fasta_competitors, sequence_from_fasta, only_identical = {}, [], {}
    fasta_identical = {}
    decoy_names = io.read_decoys(decoys) if decoys else None
    scan = {"decoys_skipped": 0, "decoys_found": set(), "long_records": []}
    if background_fasta:
        query = set().union(*group_windows.values()) if group_windows else set()
        records = scan_fasta_competitors(
            background_fasta, query, window, canonical=canonical, exclude_ids=needed,
            identical=None if keep_duplicates else configured,
            identical_out=fasta_identical, keep_ids=bg_seqs,
            keep_sequences=() if keep_duplicates else set(norm.values()),
            decoys=decoy_names or (), stats=scan)
        # a record the scan left out for a configured transcript's sequence
        records.update((rid, (seqs[conf].strip().upper(), None))
                       for rid, conf in fasta_identical.items() if rid in bg_seqs)
        for rid, (seq, hits) in records.items():
            if rid in bg_seqs:
                # one transcript in both is one column; when the FASTA holds other sequence
                # for it (another release), the column is the FASTA's, which the index has
                if seq != norm[rid]:
                    candidates[rid] = (seq, "fasta")
                    sequence_from_fasta.append(rid)
            elif hits:
                candidates[rid] = (seq, "fasta")
                fasta_competitors[rid] = len(hits)
            else:
                only_identical[rid] = seq      # kept to be told apart, not to compete
        copies = index_scope.fasta_copies(background_fasta, needed,
                                          [config.get("gene")] if config.get("gene") else ())
    else:
        copies = []

    # Salmon's default index keeps one of identical sequences, so a background sequence
    # that is a configured transcript's, or one already counted, competes with nothing:
    # through 2.4.1 that held for a FASTA record identical to a configured transcript
    # and for nothing else (issue #15).  --keep-duplicates counts every copy.
    background = {}                             # id -> the sequence of its column
    identical_to_configured = identical_to_background = identical_source = None
    if keep_duplicates:
        background = {t: seq for t, (seq, _) in candidates.items()}
    else:
        identical_to_configured, identical_to_background, identical_source = {}, {}, {}
        for rid, conf in fasta_identical.items():
            if rid not in bg_seqs:
                identical_to_configured[rid], identical_source[rid] = conf, "fasta"
        first = {}
        for t, (seq, src) in candidates.items():
            if seq in configured:
                identical_to_configured[t], identical_source[t] = configured[seq], src
            elif seq in first:
                identical_to_background[t], identical_source[t] = first[seq], src
            else:
                first[seq] = t
                background[t] = seq
        for rid, seq in only_identical.items():
            if seq in first:
                identical_to_background[rid], identical_source[rid] = first[seq], "fasta"
        identical_to_configured = dict(sorted(identical_to_configured.items()))
        identical_to_background = dict(sorted(identical_to_background.items()))
        identical_source = dict(sorted(identical_source.items()))
    for t in list(background):
        if candidates[t][1] != "fasta":
            background[t] = bg_seqs[t]          # its sequence as given, as through 2.4.1
    fasta_competitors = {t: n for t, n in fasta_competitors.items() if t in background}
    gene_columns = sorted(t for t in background if t in bg_seqs)
    bg_tracks = {t: kmer_track(s, window, canonical) for t, s in background.items()}
    bg_windows = set().union(*(set(v) for v in bg_tracks.values())) if bg_tracks else set()

    # ---- per-group report ------------------------------------------------- #
    report = {}
    positions = {}
    for g, ids in group_ids.items():
        others = set()
        for g2, ws in group_windows.items():
            if g2 != g:
                others |= ws
        shared = others | bg_windows
        uniq = group_windows[g] - shared

        best = None
        positions[g] = []
        for t in ids:
            flags = [w in uniq for w in tracks[t]]
            positions[g].extend(window_positions(flags))
            stats = coverage_stats(flags, window)
            frac = informative_fraction(
                flags, len(seqs[t]), window, read_length=read_length,
                frag_mean=frag_mean, frag_sd=frag_sd, paired=paired)
            n_reads = expected_informative_reads(
                frac, tpm, len(seqs[t]), depth=depth,
                mean_efflen=mean_efflen, frag_mean=frag_mean)
            cand = dict(stats, transcript=t, transcript_length=len(seqs[t]),
                        unique_fraction=stats["unique_length"] / max(1, len(seqs[t])),
                        informative_fraction=frac,
                        expected_informative_reads=n_reads)
            if best is None or cand["expected_informative_reads"] > best["expected_informative_reads"]:
                best = cand

        report[g] = {
            "n_transcripts": len(ids),
            "n_unique_kmers": len(uniq),
            "distinguishable": len(uniq) > 0,
            "best_transcript": best["transcript"] if best else None,
            "unique_length": best["unique_length"] if best else 0,
            "unique_fraction": best["unique_fraction"] if best else 0.0,
            "n_blocks": best["n_blocks"] if best else 0,
            "max_block_length": best["max_block_length"] if best else 0,
            "informative_fraction": best["informative_fraction"] if best else 0.0,
            "expected_informative_reads": best["expected_informative_reads"] if best else 0.0,
        }

    # ---- estimability of class totals and of their contrast --------------- #
    all_tids = needed + sorted(bg_tracks)
    all_tracks = dict(tracks)
    all_tracks.update(bg_tracks)
    A, classes = compatibility_matrix(all_tracks, all_tids)
    idx = {t: j for j, t in enumerate(all_tids)}
    # every estimand below is a functional of this one A: one decomposition serves them all
    svd = np.linalg.svd(A.T, full_matrices=False)[:2] if A.size else None

    # Expected fragments per transcript at the stated design, used as the plug-in for the
    # Poisson weights.  Flat in TPM across the gene's transcripts, which is the same
    # assumption the per-class read estimate already makes; real genes are skewed.
    all_lengths = {t: len(seqs[t]) for t in needed}
    all_lengths.update({t: len(v) for t, v in background.items()})
    theta = np.array([expected_informative_reads(
        1.0, tpm, all_lengths.get(t, 0), depth=depth,
        mean_efflen=mean_efflen, frag_mean=frag_mean) for t in all_tids])
    cov = gls_covariance(A, theta)
    class_indicator = {}

    for g, ids in group_ids.items():
        c = np.zeros(len(all_tids))
        for t in ids:
            c[idx[t]] = 1.0
        report[g].update(estimability(A, c, svd=svd))
        report[g]["label"] = "the %s class total" % g
        report[g]["coherence"] = class_coherence(tracks, ids)
        class_indicator[g] = c
        rse = gls_relative_se(cov, c, theta) if report[g]["estimable"] else float("inf")
        report[g]["gls_relative_se"] = rse
        report[g]["min_resolvable_log2fc"] = min_resolvable_log2fc(rse, n_donors)
        report[g]["beyond_linear"] = rse > LINEARISATION_LIMIT
        report[g]["verdict"], report[g]["reasons"] = _verdict(
            report[g], conditioning_tau, min_informative_reads, min_log2fc)

    c = np.zeros(len(all_tids))
    for t in group_ids[pc[0]]:
        c[idx[t]] += 1.0
    for t in group_ids[pc[1]]:
        c[idx[t]] -= 1.0
    contrast = estimability(A, c, svd=svd)
    contrast["label"] = "the class contrast"
    if contrast["estimable"]:
        raw = log_ratio_se(cov, class_indicator[pc[0]], class_indicator[pc[1]], theta)
    else:
        raw = float("inf")
    contrast["gls_relative_se"] = raw
    contrast["min_resolvable_log2fc"] = min_resolvable_log2fc(raw, n_donors)
    contrast["beyond_linear"] = raw > LINEARISATION_LIMIT
    contrast["verdict"], contrast["reasons"] = _verdict(
        contrast, conditioning_tau, None, min_log2fc)
    mean_a, mean_b, efflen_ratio = class_efflen_ratio(
        [len(seqs[t]) for t in group_ids[pc[0]]],
        [len(seqs[t]) for t in group_ids[pc[1]]], frag_mean)
    contrast["class_mean_efflen"] = {pc[0]: mean_a, pc[1]: mean_b}
    contrast["log2_efflen_ratio"] = efflen_ratio
    contrast["efflen_direction_in_band"] = (
        efflen_ratio is not None
        and abs(efflen_ratio) >= DIRECTION_MIN_ABS_LOG2_EFFLEN_RATIO)
    contrast["distinguishing_window_position"] = {
        g: position_summary(positions[g]) for g in pc[:2]}

    noise = counting_noise_floor(
        report[pc[0]]["expected_informative_reads"],
        report[pc[1]]["expected_informative_reads"],
        n_donors=n_donors)

    # The gene total -- every column of the system -- is estimable exactly when no
    # column is all zero, which happens when a transcript is shorter than ``window``.
    # That is a precondition of the system, not a judgement about the classes: it does
    # not depend on the grouping, the design or a threshold.
    gene_total = estimability(A, np.ones(len(all_tids)), svd=svd)
    gene_total["transcripts_without_windows"] = [
        t for t in all_tids if not A[:, idx[t]].any()]

    # The exit status of the CLI hangs on this, not on ``verdict``: the structural
    # verdict moves with the annotation release, while a resolvable effect size is the
    # question the experimenter actually asked.  None when no effect size was asked for.
    # An estimand past the linearisation limit is not resolved whatever its figure reads.
    primary = [report[g] for g in pc] + [contrast]
    effect_resolvable = None if min_log2fc is None else all(
        math.isfinite(e["min_resolvable_log2fc"])
        and e["min_resolvable_log2fc"] <= min_log2fc
        and not e["beyond_linear"] for e in primary)

    verdicts = [report[g]["verdict"] for g in pc] + [contrast["verdict"]]
    all_reasons = sorted({r for g in pc for r in report[g]["reasons"]}
                         | set(contrast["reasons"]))
    overall = ("not_identifiable" if "not_identifiable" in verdicts
               else "weakly_identifiable" if "weakly_identifiable" in verdicts
               else "identifiable")

    return {
        "k": k,
        "window": window,
        "canonical": canonical,
        "annotation": {"ensembl_release": config.get("ensembl_release"),
                       "fetched_release": fetched_release},
        "background": {
            "gene_id": bg_gene_id,
            "gene_transcripts": gene_columns,
            "fasta": str(background_fasta) if background_fasta else None,
            # which files: a path names a place, and the file there can change
            "fasta_sha256": io.file_sha256(background_fasta) if background_fasta else None,
            "decoys": str(decoys) if decoys else None,
            "decoys_sha256": io.file_sha256(decoys) if decoys else None,
            "decoys_listed": len(decoy_names) if decoys else None,
            "decoys_skipped": scan["decoys_skipped"] if decoys else None,
            "decoys_absent": (sorted(set(decoy_names) - scan["decoys_found"])
                              if decoys else None),
            "fasta_long_records": {str(t): n for t, n in scan["long_records"]},
            "n_background_transcripts": len(gene_columns),
            "fasta_competitors": dict(sorted(fasta_competitors.items())),
            "n_fasta_competitors": len(fasta_competitors),
            "sequence_from_fasta": sorted(sequence_from_fasta),
            "keep_duplicates": keep_duplicates,
            "identical_to_configured": identical_to_configured,
            "identical_to_background": identical_to_background,
            "identical_source": identical_source,
            "same_name_copies": copies,
        },
        "design": {"read_length": read_length, "paired": paired,
                   "frag_mean": frag_mean, "frag_sd": frag_sd, "depth": depth,
                   "mean_efflen": mean_efflen, "tpm": tpm, "n_donors": n_donors,
                   "min_log2fc": min_log2fc},
        "groups": report,
        "primary_comparison": list(pc),
        "contrast": contrast,
        "counting_noise": noise,
        "gene_total": gene_total,
        "effect_resolvable": effect_resolvable,
        "n_compatibility_classes": len(classes),
        "primary_distinguishable": all(report[g]["distinguishable"] for g in pc),
        "verdict": overall,
        "reasons": all_reasons,
    }


def run(config, k=DEFAULT_K, sequences=None, **kw):
    """Thin wrapper kept for backward compatibility; see :func:`analyze`."""
    return analyze(config, k=k, sequences=sequences, **kw)
