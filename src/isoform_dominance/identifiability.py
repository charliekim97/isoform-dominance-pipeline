"""Group-level short-read identifiability for functional isoform classes.

Short-read quantifiers apportion fragments among transcripts by solving a linear
inverse problem: the expected count of each observable fragment class is a linear
function of the transcript abundances.  A quantity is *estimable* from that system
only when its coefficient vector lies in the row space of the design matrix, and it
is estimable *usefully* only when the corresponding variance factor is small.

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

and each is checked by the textbook estimability condition together with its
generalised variance factor.  Three layers are reported, cheapest first:

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
import urllib.request

import numpy as np

ENSEMBL = "https://rest.ensembl.org"

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
def _get_text(path, timeout=30):
    req = urllib.request.Request(ENSEMBL + path, headers={"Content-Type": "text/plain"})
    return urllib.request.urlopen(req, timeout=timeout).read().decode().strip()


def fetch_cdna(transcript_id):
    """Fetch a transcript's cDNA sequence from the Ensembl REST API."""
    return _get_text("/sequence/id/%s?type=cdna" % transcript_id.split(".")[0])


def fetch_gene_transcript_ids(gene, species="homo_sapiens"):
    """Every transcript id annotated for ``gene`` (all biotypes), for use as background.

    Kept separate from :mod:`isoform_dominance.annotate`, which restricts itself to
    protein-coding transcripts: for identifiability the non-coding, retained-intron
    and NMD transcripts matter, because Salmon indexes them too.
    """
    import json

    req = urllib.request.Request(
        ENSEMBL + "/lookup/symbol/%s/%s?expand=1" % (species, gene),
        headers={"Content-Type": "application/json"},
    )
    info = json.load(urllib.request.urlopen(req, timeout=30))
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
    """
    flags = list(unique_flags)
    runs = _blocks(flags)
    # a run of r unique k-mer starts covers r + k - 1 bases
    covered = sum(r + k - 1 for _, r in runs)
    block_lengths = sorted((r + k - 1 for _, r in runs), reverse=True)
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
    ``E[count_c] = sum_t theta_t A[c, t]`` up to a shared depth factor.

    Using windows rather than full fragments is deliberately conservative: a read
    contains many windows and its compatibility set is their intersection, hence never
    larger than any single window's.  Anything this matrix says is unidentifiable is
    unidentifiable for real reads too.
    """
    sig = {}
    for tid in transcript_ids:
        for w in tracks.get(tid, []):
            sig.setdefault(w, set()).add(tid)

    counts = {}
    for owners in sig.values():
        key = frozenset(owners)
        counts[key] = counts.get(key, 0) + 1

    idx = {t: j for j, t in enumerate(transcript_ids)}
    n_windows = {t: max(1, len(tracks.get(t, []))) for t in transcript_ids}
    classes = sorted(counts, key=lambda s: (-len(s), sorted(s)))
    A = np.zeros((len(classes), len(transcript_ids)), dtype=float)
    for i, key in enumerate(classes):
        for t in key:
            # windows of t carrying this signature / total windows of t
            A[i, idx[t]] = counts[key] / n_windows[t]
    return A, [sorted(c) for c in classes]


def estimability(A, c, rcond=1e-10):
    """Is the linear functional ``c'theta`` estimable from ``E[y] = A theta``, and how well?

    **Structure.**  ``c'theta`` is estimable exactly when ``c`` lies in the row space of
    ``A``.  The residual of projecting ``c`` onto that row space, relative to ``||c||``,
    is reported as ``residual``.  This part is exact and model-free: it depends only on
    which fragment classes are observable, not on how their counts are distributed.

    **Conditioning.**  ``conditioning_factor`` is ``sqrt(c' (A'A)^+ c)``.  Read it as a
    *structural conditioning proxy*, not as the standard error of anything.  It is the
    generalised-least-squares variance factor one would get under
    ``Var(y) = sigma^2 I`` -- an assumption this package does not make and that a
    short-read quantifier does not satisfy: fragment counts are heteroskedastic, and
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
    if A.size == 0 or not np.any(c):
        return {"estimable": False, "residual": float("nan"),
                "conditioning_factor": float("inf"), "rank": 0}
    # row space of A == column space of A.T
    u, s, _ = np.linalg.svd(A.T, full_matrices=False)
    tol = max(A.shape) * (s[0] if s.size else 0.0) * rcond
    r = int((s > tol).sum())
    basis = u[:, :r]
    proj = basis @ (basis.T @ c)
    denom = np.linalg.norm(c) or 1.0
    residual = float(np.linalg.norm(c - proj) / denom)
    pinv = np.linalg.pinv(A.T @ A, rcond=rcond)
    quad = float(c @ pinv @ c)
    return {
        "estimable": residual < 1e-8,
        "residual": residual,
        "conditioning_factor": float(math.sqrt(quad)) if quad > 0 else 0.0,
        "rank": r,
    }


# --------------------------------------------------------------------------- #
# top-level analysis
# --------------------------------------------------------------------------- #
def _verdict(entry, tau, min_reads=None):
    """Grade one estimand, and say why.

    Structure first: a functional outside the row space cannot be recovered at any
    depth.  Then precision, from two independent directions -- an ill-conditioned
    contrast (large variance factor) and simple lack of informative fragments.  The
    second is why "has at least one unique k-mer" is not a usable gate: a class whose
    uniqueness is thirty junction k-mers is structurally estimable and practically
    hopeless.
    """
    reasons = []
    if not entry.get("estimable"):
        reasons.append("contrast outside the row space of the compatibility system")
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
            min_informative_reads=DEFAULT_MIN_INFORMATIVE_READS):
    """Assess whether the configured isoform classes are measurable by short reads.

    Parameters
    ----------
    config
        The pipeline config dict; ``groups`` and ``primary_comparison`` are used.
    k, window
        k-mer length, and the window length used to build the compatibility system
        (defaults to ``k``; set it to the read length for a sharper, still
        conservative, system).
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
        configured groups only.
    read_length, frag_mean, frag_sd, paired, depth, mean_efflen, tpm, n_donors
        The sequencing design the report should be conditioned on.
    conditioning_tau, min_informative_reads
        Thresholds separating ``identifiable`` from ``weakly_identifiable``.

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
    if background_gene_transcripts and config.get("gene") and not bg_seqs:
        try:
            all_ids = fetch_gene_transcript_ids(
                config["gene"], species or config.get("species", "homo_sapiens"))
            for tid in all_ids:
                if tid not in needed and tid not in bg_seqs:
                    bg_seqs[tid] = fetch_cdna(tid)
        except Exception:      # offline, or symbol not found: fall back quietly
            bg_seqs = dict(bg_seqs)
    for tid in needed:
        if tid not in seqs:
            seqs[tid] = fetch_cdna(tid)
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
        report[g]["verdict"], report[g]["reasons"] = _verdict(
            report[g], conditioning_tau, min_informative_reads)

    c = np.zeros(len(all_tids))
    for t in group_ids[pc[0]]:
        c[idx[t]] += 1.0
    for t in group_ids[pc[1]]:
        c[idx[t]] -= 1.0
    contrast = estimability(A, c)
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
