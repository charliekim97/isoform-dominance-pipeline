"""The v2.2 identifiability layers: background, canonical k-mers, read model, estimability.

Every test is offline and deterministic: sequences are supplied explicitly and
``background_gene_transcripts`` therefore never reaches for the network.
"""
import gzip
import random

import numpy as np
import pytest

from isoform_dominance import identifiability as I


def _seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


SHARED = _seq(1200, 1)
EXON_A = _seq(700, 2)
EXON_B = _seq(900, 3)


def _cfg():
    return {"gene": "TEST",
            "groups": {"short": ["A1"], "long": ["B1"]},
            "primary_comparison": ["short", "long"]}


def _seqs():
    return {"A1": SHARED + EXON_A, "B1": SHARED + EXON_B}


# --------------------------------------------------------------------------- #
# canonical k-mers
# --------------------------------------------------------------------------- #
def test_revcomp_and_canonical_kmer():
    assert I.revcomp("ACGTN") == "NACGT"
    assert I.canonical_kmer("AAAA") == "AAAA"          # rc TTTT, AAAA is smaller
    assert I.canonical_kmer("TTTT") == "AAAA"          # folds onto its rc


def test_canonical_folding_removes_orientation_only_uniqueness():
    """B carries the reverse complement of A's "unique" stretch.

    Strand-aware counting calls the whole stretch unique to each class; a canonical
    index -- which is what Salmon builds, and what an unstranded library requires --
    cannot tell the two orientations apart.  Only the k-mers straddling the
    shared/tail junction, which differ between the two transcripts as sequence,
    survive canonicalisation.
    """
    tail = _seq(120, 7)
    seqs = {"A1": SHARED + tail, "B1": SHARED + I.revcomp(tail)}
    cfg = {"groups": {"a": ["A1"], "b": ["B1"]}, "primary_comparison": ["a", "b"]}
    strand_aware = I.analyze(cfg, k=31, sequences=seqs, canonical=False,
                             background_gene_transcripts=False)
    canonical = I.analyze(cfg, k=31, sequences=seqs, canonical=True,
                          background_gene_transcripts=False)
    assert strand_aware["groups"]["a"]["n_unique_kmers"] >= 100
    assert canonical["groups"]["a"]["n_unique_kmers"] < 31        # junction only
    assert (canonical["groups"]["a"]["n_unique_kmers"]
            < strand_aware["groups"]["a"]["n_unique_kmers"] / 4)


# --------------------------------------------------------------------------- #
# background
# --------------------------------------------------------------------------- #
def test_background_sequences_remove_false_uniqueness():
    """A k-mer unique among the configured groups is not unique in the index.

    Up to v2.1.1 uniqueness was judged only against the other configured groups, so an
    unlisted transcript of the same gene (here a retained-intron form carrying the
    whole short-class terminal exon) could not lower the verdict.
    """
    naive = I.analyze(_cfg(), k=31, sequences=_seqs(),
                      background_gene_transcripts=False)
    assert naive["groups"]["short"]["n_unique_kmers"] == len(EXON_A)

    withbg = I.analyze(_cfg(), k=31, sequences=_seqs(),
                       background_sequences={"RI1": _seq(300, 9) + EXON_A},
                       background_gene_transcripts=False)
    # only the k-mers spanning the shared/exon junction remain unique
    assert withbg["groups"]["short"]["n_unique_kmers"] < 31
    assert withbg["background"]["n_background_transcripts"] == 1
    # structurally still estimable, but now starved of informative fragments --
    # which is exactly the case a bare "has a unique k-mer" gate waves through
    assert withbg["groups"]["short"]["estimable"] is True
    assert withbg["verdict"] == "weakly_identifiable"
    assert any("informative reads" in r for r in withbg["reasons"])


@pytest.mark.parametrize("gz", [False, True])
def test_scan_background_fasta_streams_plain_and_gzip(tmp_path, gz):
    path = tmp_path / ("bg.fa.gz" if gz else "bg.fa")
    body = ">ENST00000000001.4|extra|fields\n%s\n>ENST00000000002\n%s\n" % (
        EXON_A, _seq(200, 11))
    if gz:
        path.write_bytes(gzip.compress(body.encode()))
    else:
        path.write_text(body)

    query = I.kmers(EXON_A, 31)
    hits = I.scan_background_fasta(path, query, 31)
    assert hits == query                                   # every EXON_A k-mer is present

    # a record named in exclude_ids is not allowed to be its own background
    hits = I.scan_background_fasta(path, query, 31,
                                   exclude_ids=["ENST00000000001"])
    assert hits == set()


def test_background_fasta_lowers_the_verdict(tmp_path):
    path = tmp_path / "bg.fa"
    path.write_text(">DECOY\n%s\n" % (SHARED + EXON_A))
    res = I.analyze(_cfg(), k=31, sequences=_seqs(), background_fasta=str(path),
                    background_gene_transcripts=False)
    assert res["groups"]["short"]["n_unique_kmers"] == 0
    assert res["background"]["fasta"].endswith("bg.fa")


def test_auto_background_stays_offline_when_sequences_are_supplied():
    # "auto" must not touch the network when the caller already supplied every
    # sequence -- otherwise an offline run blocks on a 30 s timeout per transcript.
    res = I.analyze(_cfg(), k=31, sequences=_seqs())      # background_gene_transcripts="auto"
    assert res["background"]["n_background_transcripts"] == 0


# --------------------------------------------------------------------------- #
# coverage / blocks
# --------------------------------------------------------------------------- #
def test_coverage_stats_counts_bases_not_kmers():
    # 5 consecutive unique k-mer starts at k=4 cover 5 + 4 - 1 = 8 bases
    flags = [False, True, True, True, True, True, False, False]
    st = I.coverage_stats(flags, 4)
    assert st == {"unique_length": 8, "n_blocks": 1, "max_block_length": 8,
                  "block_lengths": [8]}


def test_coverage_stats_splits_blocks():
    flags = [True, True, False, False, True]
    st = I.coverage_stats(flags, 3)
    assert st["n_blocks"] == 2
    assert st["max_block_length"] == 4          # 2 starts + 3 - 1
    assert st["unique_length"] == 4 + 3


def test_unique_region_is_contiguous_for_an_alternative_terminal_exon():
    res = I.analyze(_cfg(), k=31, sequences=_seqs(),
                    background_gene_transcripts=False)
    short = res["groups"]["short"]
    assert short["n_blocks"] == 1
    # the unique run starts k-1 bases before the junction, so the covered span is
    # the alternative exon plus the k-1 bases of shared sequence it is read against
    assert short["unique_length"] == len(EXON_A) + 31 - 1
    assert 0.0 < short["unique_fraction"] < 1.0


# --------------------------------------------------------------------------- #
# read / fragment model
# --------------------------------------------------------------------------- #
def test_informative_fraction_ignores_unique_sequence_only_reachable_mid_fragment():
    """Reads, not fragments, carry the evidence.

    A unique island buried in the middle of a long fragment is never sequenced when
    the reads are short, so it must not count towards informativeness.
    """
    n_starts = 2000
    flags = [False] * n_starts
    flags[950:1050] = [True] * 100               # island near the transcript's middle
    reachable = I.informative_fraction(flags, 2100, 31, read_length=100,
                                       frag_mean=200, frag_sd=1, paired=True)
    # with a 2000 nt fragment the island sits far from both sequenced ends
    unreachable = I.informative_fraction(flags, 2100, 31, read_length=100,
                                         frag_mean=2000, frag_sd=1, paired=True)
    assert reachable > 0.05
    assert unreachable < reachable / 5


def test_informative_fraction_is_zero_without_unique_kmers():
    assert I.informative_fraction([False] * 500, 530, 31) == 0.0


def test_paired_end_sees_at_least_as_much_as_single_end():
    flags = [False] * 1000
    flags[900:1000] = [True] * 100               # unique stretch at the 3' end
    single = I.informative_fraction(flags, 1030, 31, paired=False)
    paired = I.informative_fraction(flags, 1030, 31, paired=True)
    assert paired >= single


def test_expected_informative_reads_scales_with_depth_and_tpm():
    base = I.expected_informative_reads(0.5, tpm=10, transcript_length=2000)
    assert I.expected_informative_reads(0.5, tpm=20, transcript_length=2000) == pytest.approx(2 * base)
    assert I.expected_informative_reads(0.5, tpm=10, transcript_length=2000,
                                        depth=60_000_000) == pytest.approx(2 * base)


def test_counting_noise_floor_is_nan_without_informative_reads():
    assert np.isnan(I.counting_noise_floor(0, 100)["log2_ratio_se"])
    got = I.counting_noise_floor(100, 100, n_donors=4)
    assert got["log2_ratio_se"] > 0
    # averaging over donors tightens the bound by sqrt(n)
    assert got["min_resolvable_log2fc"] == pytest.approx(
        1.96 * got["log2_ratio_se"] / 2.0)


# --------------------------------------------------------------------------- #
# estimability
# --------------------------------------------------------------------------- #
def test_nested_class_without_unique_kmers_is_still_estimable():
    """Zero unique k-mers does not imply unmeasurable -- the v2.1 gate's false negative.

    Here the short class is a strict prefix of the long one and owns no unique k-mer
    at all, so the ``n_unique_kmers > 0`` test refuses the comparison.  The class total
    is nonetheless an estimable function of the system: the containing class is pinned
    by its own unique region, so the nested class's indicator lies in the row space.

    Note what is *not* being claimed -- that an EM quantifier "works it out".  EM
    returns numbers for a non-identifiable model too; the likelihood is simply flat
    along the unidentified directions.  Estimability is a property of the design
    matrix, and that is what is asserted here.
    """
    cfg = {"groups": {"sub": ["S1"], "sup": ["P1"]},
           "primary_comparison": ["sub", "sup"]}
    res = I.analyze(cfg, k=31, sequences={"S1": SHARED, "P1": SHARED + EXON_B},
                    background_gene_transcripts=False)
    assert res["groups"]["sub"]["n_unique_kmers"] == 0
    assert res["groups"]["sub"]["distinguishable"] is False      # the old verdict
    assert res["groups"]["sub"]["estimable"] is True             # the correct one
    assert res["contrast"]["estimable"] is True


def test_identical_transcripts_break_each_class_but_not_their_sum():
    """Individually unidentifiable, jointly identifiable -- the group-level point.

    Two transcripts with identical sequence produce a single compatibility class, so
    neither abundance nor their difference is estimable, while the total is.  A tool
    that grouped these two transcripts together would be measuring something real;
    one that contrasts them is measuring noise.
    """
    dup = SHARED + EXON_A
    cfg = {"groups": {"g1": ["D1"], "g2": ["D2"]},
           "primary_comparison": ["g1", "g2"]}
    res = I.analyze(cfg, k=31, sequences={"D1": dup, "D2": dup},
                    background_gene_transcripts=False)
    assert res["groups"]["g1"]["estimable"] is False
    assert res["groups"]["g2"]["estimable"] is False
    assert res["contrast"]["estimable"] is False
    assert res["verdict"] == "not_identifiable"

    tracks = {"D1": I.kmer_track(dup, 31), "D2": I.kmer_track(dup, 31)}
    A, _ = I.compatibility_matrix(tracks, ["D1", "D2"])
    assert I.estimability(A, np.array([1.0, 1.0]))["estimable"] is True
    assert I.estimability(A, np.array([1.0, -1.0]))["estimable"] is False


def test_alternative_terminal_exons_are_estimable_with_a_small_conditioning_factor():
    res = I.analyze(_cfg(), k=31, sequences=_seqs(),
                    background_gene_transcripts=False)
    for g in ("short", "long"):
        assert res["groups"][g]["estimable"] is True
        assert res["groups"][g]["conditioning_factor"] < I.DEFAULT_CONDITIONING_TAU
        assert res["groups"][g]["verdict"] == "identifiable"
    assert res["contrast"]["estimable"] is True
    assert res["verdict"] == "identifiable"


def test_conditioning_factor_threshold_downgrades_to_weakly_identifiable():
    res = I.analyze(_cfg(), k=31, sequences=_seqs(),
                    background_gene_transcripts=False,
                    conditioning_tau=1e-6)
    assert res["verdict"] == "weakly_identifiable"
    assert res["groups"]["short"]["estimable"] is True
    assert any("conditioning factor" in r for r in res["reasons"])


def test_the_zero_functional_is_estimable():
    """0'theta = 0 lies in every row space and is estimated by the constant 0.

    The previous answer here was `estimable=False, conditioning=inf, rank=0`, which is
    wrong on all three counts and was reachable from a config whose two groups name the
    same transcripts: the tool refused the self-comparison, correctly, for a reason that
    was not true.
    """
    A = np.eye(2)
    got = I.estimability(A, np.zeros(2))
    assert got["estimable"] is True
    assert got["residual"] == 0.0
    assert got["conditioning_factor"] == 0.0
    assert got["rank"] == 2

    # with no observable classes the row space is {0}: only c = 0 survives
    empty = np.zeros((0, 2))
    assert I.estimability(empty, np.zeros(2))["estimable"] is True
    assert I.estimability(empty, np.array([1.0, -1.0]))["estimable"] is False


def test_conditioning_is_truncated_on_the_same_tolerance_as_the_rank():
    """A contrast that is declared estimable must get a conditioning factor to match.

    Regression for a tolerance mismatch: the rank test ran on sigma(A) while the
    conditioning ran through `pinv(A'A)`, whose rcond is relative to sigma(A)**2. Every
    direction in the five-decade band between the two thresholds was declared estimable
    and had its conditioning direction projected away at the same time, so the reported
    factor collapsed to 0.0 -- the best possible score -- for the worst-conditioned
    contrasts there are. That inverts the reason the quantity exists.
    """
    for eps in (1e-3, 1e-5, 1e-7, 1e-9):
        got = I.estimability(np.diag([1.0, eps]), np.array([0.0, 1.0]))
        assert got["estimable"] is True
        assert got["rank"] == 2
        # exact answer for a diagonal design is 1/eps
        assert got["conditioning_factor"] == pytest.approx(1.0 / eps, rel=1e-6)

    # and a contrast outside the row space has no finite factor at all
    bad = I.estimability(np.array([[1.0, 1.0]]), np.array([1.0, -1.0]))
    assert bad["estimable"] is False
    assert np.isinf(bad["conditioning_factor"])


def test_compatibility_matrix_collapses_indistinguishable_windows():
    A, classes = I.compatibility_matrix(
        {"t1": ["aa", "bb"], "t2": ["bb", "cc"]}, ["t1", "t2"])
    # three signatures: {t1}, {t1,t2}, {t2}
    assert len(classes) == 3
    assert A.shape == (3, 2)
    # rows are per-transcript probabilities, so each column sums to 1
    assert A.sum(axis=0) == pytest.approx([1.0, 1.0])


def test_compatibility_matrix_counts_positions_not_distinct_windows():
    """A repeated window is drawn twice as often, so it must be counted twice.

    Counting distinct window *sequences* per signature while dividing by the number of
    window *positions* breaks the column-sum invariant and under-weights precisely the
    shared classes that absorb the most fragments -- and repeated windows are the rule
    in cDNA (A-rich 3' ends, tandem repeats, Alu in long UTRs), not an edge case. The
    previous fixture could not see it because all of its windows were distinct.
    """
    A, classes = I.compatibility_matrix(
        {"t1": ["aa", "bb", "aa"], "t2": ["bb", "cc"]}, ["t1", "t2"])
    assert A.sum(axis=0) == pytest.approx([1.0, 1.0])
    own = classes.index(["t1"])
    shared = classes.index(["t1", "t2"])
    # two of t1's three window positions are its own; one is shared
    assert A[own, 0] == pytest.approx(2.0 / 3.0)
    assert A[shared, 0] == pytest.approx(1.0 / 3.0)

    # a real repeat: an A-rich 3' end repeats one 31-mer many times over
    tracks = {"x": I.kmer_track("ACGT" * 40 + "A" * 80, 31, True),
              "y": I.kmer_track("ACGT" * 40 + "TTTT" * 20, 31, True)}
    B, _ = I.compatibility_matrix(tracks, ["x", "y"])
    assert B.sum(axis=0) == pytest.approx([1.0, 1.0])


# --------------------------------------------------------------------------- #
# report shape
# --------------------------------------------------------------------------- #
def test_report_keeps_the_v21_public_keys():
    res = I.analyze(_cfg(), k=31, sequences=_seqs(),
                    background_gene_transcripts=False)
    for key in ("k", "groups", "primary_comparison", "primary_distinguishable"):
        assert key in res
    for key in ("n_unique_kmers", "n_transcripts", "distinguishable"):
        assert key in res["groups"]["short"]


def test_window_may_exceed_k_for_a_sharper_system():
    res = I.analyze(_cfg(), k=31, sequences=_seqs(), window=100,
                    background_gene_transcripts=False)
    assert res["window"] == 100
    assert res["groups"]["short"]["estimable"] is True


def test_empty_groups_rejected():
    with pytest.raises(ValueError, match="empty"):
        I.analyze({"groups": {}, "primary_comparison": []}, sequences={})


# --------------------------------------------------------------------------- #
# The README and paper claim that transcript-level identifiability is neither
# necessary nor sufficient for a class contrast to be measurable.  The
# "not necessary" half is covered by
# test_identical_transcripts_break_each_class_but_not_their_sum above.  These two
# pin the "not sufficient" half, and pin it precisely, because the obvious
# statement of it is false: full column rank does imply that every contrast is
# estimable.  What it does not imply is that any of them is recoverable.
def test_full_column_rank_makes_every_contrast_estimable():
    """Individually identifiable transcripts cannot yield a non-estimable contrast.

    If every unit vector lies in the row space then the row space is the whole
    space, so no functional is left out.  Stating otherwise -- that identifiable
    transcripts can define a contrast that is not estimable -- is wrong, and a
    reviewer refutes it in one line.  The insufficiency is about conditioning,
    not about the row-space condition; the next test carries that.
    """
    A = np.array([[1.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0],
                  [0.0, 0.0, 1.0],
                  [1.0, 1.0, 1.0]])
    assert I.estimability(A, np.array([1.0, 1.0, 1.0]))["rank"] == 3
    for t in range(3):
        assert I.estimability(A, np.eye(3)[t])["estimable"] is True
    rng = np.random.default_rng(0)
    for _ in range(50):
        c = rng.normal(size=3)
        assert I.estimability(A, c)["estimable"] is True


def test_full_column_rank_does_not_bound_the_conditioning_factor():
    """Every transcript identifiable, the contrast estimable, and still unusable.

    Here the class direction is observed only through a difference of nearly
    identical rows, so the contrast is estimable in exact arithmetic and its
    conditioning factor is three orders of magnitude worse than a well-posed
    design.  This is why the rank verdict is reported with the conditioning
    factor rather than on its own.
    """
    eps = 1e-3
    A = np.array([[1.0, 1.0, 1.0],
                  [1.0, 1.0, 1.0 + eps],
                  [1.0, 1.0 - eps, 1.0]])
    c = np.array([1.0, 1.0, -1.0])
    got = I.estimability(A, c)
    assert got["rank"] == 3
    for t in range(3):
        assert I.estimability(A, np.eye(3)[t])["estimable"] is True
    assert got["estimable"] is True
    assert got["conditioning_factor"] > 1e3

    well_posed = np.eye(3)
    assert I.estimability(well_posed, c)["conditioning_factor"] < 2.0
