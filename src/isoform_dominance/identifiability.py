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
from urllib.error import HTTPError

import numpy as np

from . import ensembl
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


def fetch_gene_transcript_ids(gene, species="homo_sapiens", **retry):
    """Every transcript id annotated for ``gene`` (all biotypes), for use as background.

    Kept separate from :mod:`isoform_dominance.annotate`, which restricts itself to
    protein-coding transcripts: for identifiability the non-coding, retained-intron
    and NMD transcripts matter, because Salmon indexes them too.
    """
    info = ensembl.get_json("/lookup/symbol/%s/%s?expand=1" % (species, gene), **retry)
    return [t["id"].split(".")[0] for t in info.get("Transcript", [])]


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
    seq = seq.upper()
    if k <= 0 or len(seq) < k:
        return set()
    if canonical:
        return {canonical_kmer(seq[i:i + k]) for i in range(len(seq) - k + 1)}
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def kmer_track(seq, k, canonical=True):
    """Ordered list of the k-mers of ``seq``, one per start position.

    Unlike :func:`kmers` this preserves position, which is what the coverage,
    block-structure and read-model calculations need.
    """
    seq = seq.upper()
    if k <= 0 or len(seq) < k:
        return []
    if canonical:
        return [canonical_kmer(seq[i:i + k]) for i in range(len(seq) - k + 1)]
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def scan_background_fasta(path, query_kmers, k, canonical=True, exclude_ids=()):
    """Return the subset of ``query_kmers`` that also occurs in a background FASTA.

    Streams the file and never materialises the background's own k-mer set, so a
    whole-transcriptome FASTA (GENCODE, or the FASTA a Salmon index was built from)
    can be used as background in memory proportional to the *query*, not the file.
    Handles plain or gzipped input; ``exclude_ids`` drops records whose first
    ``|``- or whitespace-delimited field matches (version suffix ignored), which is
    how the transcripts under test are kept out of their own background.
    """
    import gzip

    query = set(query_kmers)
    if not query:
        return set()
    exclude = {str(i).split(".")[0] for i in exclude_ids}
    seen = set()
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open

    def _consume(chunks, tid):
        if not chunks or tid in exclude:
            return
        seq = "".join(chunks).upper()
        for i in range(len(seq) - k + 1):
            km = seq[i:i + k]
            if canonical:
                km = canonical_kmer(km)
            if km in query:
                seen.add(km)

    chunks, tid = [], None
    with opener(path, "rt") as fh:
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
    efflen = max(1.0, transcript_length - frag_mean + 1.0)
    reads = depth * (tpm * efflen) / (1e6 * max(1.0, mean_efflen))
    return float(reads * informative_frac)


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
    sig = {}
    for tid in transcript_ids:
        for w in tracks.get(tid, []):
            sig.setdefault(w, set()).add(tid)

    # positions carrying each signature, per transcript -- see the docstring on why
    # this is not a count of distinct window sequences
    pos_counts = {}
    class_keys = set()
    for tid in transcript_ids:
        for w in tracks.get(tid, []):
            key = frozenset(sig[w])
            class_keys.add(key)
            pos_counts[(key, tid)] = pos_counts.get((key, tid), 0) + 1

    idx = {t: j for j, t in enumerate(transcript_ids)}
    n_windows = {t: max(1, len(tracks.get(t, []))) for t in transcript_ids}
    classes = sorted(class_keys, key=lambda s: (-len(s), sorted(s)))
    A = np.zeros((len(classes), len(transcript_ids)), dtype=float)
    for i, key in enumerate(classes):
        for t in key:
            A[i, idx[t]] = pos_counts.get((key, t), 0) / n_windows[t]
    return A, [sorted(c) for c in classes]


def estimability(A, c, rcond=1e-10):
    """Is the linear functional ``c'theta`` estimable from ``E[y] = A theta``, and how well?

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
    u, s, _ = np.linalg.svd(A.T, full_matrices=False)
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
def _verdict(entry, tau, min_reads=None):
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
    if entry.get("conditioning_factor", 0.0) > tau:
        reasons.append("conditioning factor %.1f exceeds tau=%.1f"
                       % (entry["conditioning_factor"], tau))
    n = entry.get("expected_informative_reads")
    if min_reads is not None and n is not None and n < min_reads:
        reasons.append("expected informative reads %.0f below %.0f" % (n, min_reads))
    return ("weakly_identifiable" if reasons else "identifiable"), reasons


def analyze(config, k=DEFAULT_K, sequences=None, *,
            canonical=True, window=None,
            background_sequences=None, background_fasta=None,
            background_gene_transcripts="auto", species=None,
            read_length=DEFAULT_READ_LENGTH, frag_mean=DEFAULT_FRAG_MEAN,
            frag_sd=DEFAULT_FRAG_SD, paired=True, depth=DEFAULT_DEPTH,
            mean_efflen=DEFAULT_MEAN_EFFLEN, tpm=DEFAULT_TPM, n_donors=1,
            conditioning_tau=DEFAULT_CONDITIONING_TAU,
            min_informative_reads=DEFAULT_MIN_INFORMATIVE_READS,
            retries=DEFAULT_RETRIES, retry_wait=DEFAULT_RETRY_WAIT):
    """Assess whether the configured isoform classes are measurable by short reads.

    Parameters
    ----------
    config
        The pipeline config dict; ``groups`` and ``primary_comparison`` are used.
    k, window
        k-mer length, and the window length used to build the compatibility system
        (defaults to ``k``; a different window gives a different system, not a
        uniformly sharper one -- the rank is not monotone in it, see
        :func:`compatibility_matrix`)
    sequences
        ``{transcript_id: cdna}``.  Anything missing is fetched from Ensembl.
    background_sequences, background_fasta, background_gene_transcripts
        What uniqueness is judged against, beyond the other configured groups.  Pass a
        FASTA -- ideally the one the Salmon index was built from -- for the honest
        whole-index answer.  ``background_gene_transcripts`` defaults to ``"auto"``:
        the gene's remaining transcripts are fetched and used when the caller is
        already relying on Ensembl for sequence, and skipped when sequences were
        supplied offline (so an offline call never blocks on the network).  ``True``
        forces the fetch, ``False`` restores the pre-v2.2 behaviour of comparing the
        configured groups only.  A gene symbol Ensembl does not know (HTTP 400/404)
        leaves the gene background empty; any other failure to fetch it is raised,
        because a background that is only partly fetched silently gives a different
        answer.
    read_length, frag_mean, frag_sd, paired, depth, mean_efflen, tpm, n_donors
        The sequencing design the report should be conditioned on.
    conditioning_tau, min_informative_reads
        Thresholds separating ``identifiable`` from ``weakly_identifiable``.
    retries, retry_wait
        Each Ensembl request is retried up to ``retries`` times after the first
        attempt, waiting ``retry_wait`` seconds and doubling each time (an HTTP 429
        waits for its ``Retry-After``); defaults 5 and 1.0.  cDNA is fetched 50
        transcripts per request.  See :mod:`isoform_dominance.ensembl`.

    Returns
    -------
    dict
        ``groups`` (per-class sequence, read-model and estimability numbers),
        ``primary_comparison``, ``contrast`` (estimability of ``s_A - s_B``),
        ``verdict`` with its ``reasons``, and -- for callers written against v2.1 --
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

    pc = config.get("primary_comparison", list(groups)[:2])
    if len(pc) < 2:
        raise ValueError(
            "primary_comparison must name two isoform groups; got %r" % (pc,))
    missing = [g for g in pc if g not in groups]
    if missing:
        raise ValueError(
            "primary_comparison names group(s) not in config['groups']: %r" % (missing,))

    group_ids = {g: [t.split(".")[0] for t in ids] for g, ids in groups.items()}
    needed = [t for ids in group_ids.values() for t in ids]
    seqs = {tid.split(".")[0]: s for tid, s in (sequences or {}).items()}

    # ---- background transcripts of the same gene -------------------------- #
    bg_seqs = {t.split(".")[0]: s for t, s in (background_sequences or {}).items()}
    if background_gene_transcripts == "auto":
        # only reach for the network when we are already going there for sequence
        background_gene_transcripts = any(t not in seqs for t in needed)
    net = {"retries": retries, "retry_wait": retry_wait}
    if background_gene_transcripts and config.get("gene") and not bg_seqs:
        try:
            all_ids = fetch_gene_transcript_ids(
                config["gene"], species or config.get("species", "homo_sapiens"), **net)
        except HTTPError as e:
            if e.code not in (400, 404):      # 400 is Ensembl's "no such symbol"
                raise
            all_ids = []                      # a symbol Ensembl does not know
        # all or nothing: this used to swallow any error part-way through and carry
        # on with whatever had arrived, which is a different answer, silently
        bg_seqs.update(ensembl.fetch_cdna_batch(
            [t for t in all_ids if t not in needed], **net))
    missing = [t for t in needed if t not in seqs]
    if missing:
        got = ensembl.fetch_cdna_batch(missing, **net)
        absent = [t for t in missing if t not in got]
        if absent:
            raise ValueError("Ensembl returned no cDNA for %s" % ", ".join(absent))
        seqs.update(got)
    bg_seqs = {t: s for t, s in bg_seqs.items() if t not in needed}

    # ---- window tracks ---------------------------------------------------- #
    tracks = {t: kmer_track(seqs[t], window, canonical) for t in needed}
    bg_tracks = {t: kmer_track(s, window, canonical) for t, s in bg_seqs.items()}

    group_windows = {g: set().union(*(set(tracks[t]) for t in ids)) if ids else set()
                     for g, ids in group_ids.items()}
    bg_windows = set().union(*(set(v) for v in bg_tracks.values())) if bg_tracks else set()

    # ---- an external FASTA background, streamed --------------------------- #
    fasta_hits = set()
    if background_fasta:
        query = set().union(*group_windows.values()) if group_windows else set()
        fasta_hits = scan_background_fasta(
            background_fasta, query, window, canonical=canonical, exclude_ids=needed)

    # ---- per-group report ------------------------------------------------- #
    report = {}
    for g, ids in group_ids.items():
        others = set()
        for g2, ws in group_windows.items():
            if g2 != g:
                others |= ws
        shared = others | bg_windows | fasta_hits
        uniq = group_windows[g] - shared

        best = None
        for t in ids:
            flags = [w in uniq for w in tracks[t]]
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

    for g, ids in group_ids.items():
        c = np.zeros(len(all_tids))
        for t in ids:
            c[idx[t]] = 1.0
        report[g].update(estimability(A, c))
        report[g]["label"] = "the %s class total" % g
        report[g]["verdict"], report[g]["reasons"] = _verdict(
            report[g], conditioning_tau, min_informative_reads)

    c = np.zeros(len(all_tids))
    for t in group_ids[pc[0]]:
        c[idx[t]] += 1.0
    for t in group_ids[pc[1]]:
        c[idx[t]] -= 1.0
    contrast = estimability(A, c)
    contrast["label"] = "the class contrast"
    contrast["verdict"], contrast["reasons"] = _verdict(contrast, conditioning_tau)

    noise = counting_noise_floor(
        report[pc[0]]["expected_informative_reads"],
        report[pc[1]]["expected_informative_reads"],
        n_donors=n_donors)

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
        "background": {
            "gene_transcripts": sorted(bg_tracks),
            "fasta": str(background_fasta) if background_fasta else None,
            "n_background_transcripts": len(bg_tracks),
        },
        "design": {"read_length": read_length, "paired": paired,
                   "frag_mean": frag_mean, "frag_sd": frag_sd, "depth": depth,
                   "mean_efflen": mean_efflen, "tpm": tpm, "n_donors": n_donors},
        "groups": report,
        "primary_comparison": list(pc),
        "contrast": contrast,
        "counting_noise": noise,
        "n_compatibility_classes": len(classes),
        "primary_distinguishable": all(report[g]["distinguishable"] for g in pc),
        "verdict": overall,
        "reasons": all_reasons,
    }


def run(config, k=DEFAULT_K, sequences=None, **kw):
    """Thin wrapper kept for backward compatibility; see :func:`analyze`."""
    return analyze(config, k=k, sequences=sequences, **kw)
