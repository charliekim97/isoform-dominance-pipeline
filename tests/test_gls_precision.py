"""The Poisson-GLS precision figures, and the properties that make them trustworthy.

Three claims are load-bearing and each is checked against something outside the formula:

* the GLS estimate is the BLUE, so it cannot be beaten by an estimator built from the
  compatibility rows that only one class contributes to;
* the delta-method figure for ``log2(A/B)`` matches simulation while the linearisation
  holds, and is flagged when it does not;
* class coherence separates a group whose members resemble each other from one that
  merely shares a terminal acceptor.
"""
import math
import random

import numpy as np
import pytest

from isoform_dominance import identifiability as I

W = 31


def _rng_seqs(spec, seed=0):
    rng = random.Random(seed)

    def draw(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    core_a, core_b = draw(700), draw(700)
    out = {}
    for name, (which, private) in spec.items():
        out[name] = (core_a if which == "A" else core_b) + draw(private)
    return out


def _system(seqs, groups, window=W):
    ids = [t for v in groups.values() for t in v]
    tracks = {t: I.kmer_track(seqs[t], window, True) for t in ids}
    out = I.compatibility_matrix(tracks, ids)
    A = out[0] if isinstance(out, tuple) else out
    ind = {g: np.array([1.0 if t in v else 0.0 for t in ids]) for g, v in groups.items()}
    return A, ind, ids, tracks


def test_gls_beats_any_estimator_built_from_one_class_own_rows():
    """BLUE: using every fragment cannot be worse than using the unambiguous ones.

    The comparison estimator takes the compatibility rows no other class contributes to
    and rescales by the share of the class's mass that lands there -- unbiased, and the
    natural formalisation of "count the reads that can only have come from this class".
    """
    seqs = _rng_seqs({"A0": ("A", 400), "A1": ("A", 300), "A2": ("A", 350),
                      "B0": ("B", 400), "B1": ("B", 300)})
    groups = {"A": ["A0", "A1", "A2"], "B": ["B0", "B1"]}
    A, ind, ids, _ = _system(seqs, groups)
    theta = np.full(len(ids), 2000.0)
    mean = A @ theta
    cov = I.gls_covariance(A, theta)
    assert cov is not None

    for g, c in ind.items():
        est = float(c @ theta)
        own = A @ c
        rest = A @ (np.ones(len(ids)) - c)
        only = (own > 1e-12) & (rest < 1e-12)
        assert only.any(), "this fixture must give each class rows of its own"
        share = float(own[only].sum()) / est
        var_own_rows = float(mean[only].sum()) / (share * share)
        var_gls = float(c @ cov @ c)
        assert var_gls <= var_own_rows * (1 + 1e-9), (
            "GLS variance %g exceeded the own-rows estimator's %g for class %s"
            % (var_gls, var_own_rows, g))


def test_the_log_ratio_se_matches_simulation_while_the_linearisation_holds():
    seqs = _rng_seqs({"A0": ("A", 400), "A1": ("A", 300),
                      "B0": ("B", 400), "B1": ("B", 300)}, seed=3)
    groups = {"A": ["A0", "A1"], "B": ["B0", "B1"]}
    A, ind, ids, _ = _system(seqs, groups)
    theta = np.full(len(ids), 20000.0)
    mean = A @ theta
    cov = I.gls_covariance(A, theta)
    se = I.log_ratio_se(cov, ind["A"], ind["B"], theta)
    assert se < I.LINEARISATION_LIMIT, "fixture must sit inside the linear regime"

    keep = mean > 1e-12
    aw = A[keep] / np.sqrt(mean[keep])[:, None]
    pinv = np.linalg.pinv(aw.T @ aw, rcond=1e-10)
    rng = np.random.default_rng(11)
    draws = []
    for _ in range(3000):
        y = rng.poisson(mean)
        fit = pinv @ (aw.T @ (y[keep] / np.sqrt(mean[keep])))
        a, b = float(ind["A"] @ fit), float(ind["B"] @ fit)
        if a > 0 and b > 0:
            draws.append(math.log(a / b))
    empirical = float(np.std(draws, ddof=1))
    assert 0.8 < empirical / se < 1.25, (
        "delta-method SE %g against simulated %g" % (se, empirical))


def test_a_figure_past_the_linearisation_limit_is_flagged():
    """Starve the design until the first-order figure stops being a value.

    Simulation puts the delta-method SE within a few percent of the truth below
    :data:`I.LINEARISATION_LIMIT` and up to twofold away above it, so the flag is what
    separates a number to quote from one that only means "out of reach".
    """
    seqs = _rng_seqs({"A0": ("A", 400), "A1": ("A", 300),
                      "B0": ("B", 400), "B1": ("B", 300)}, seed=5)
    cfg = {"groups": {"A": ["A0", "A1"], "B": ["B0", "B1"]},
           "primary_comparison": ["A", "B"]}

    rich = I.analyze(cfg, sequences=seqs, background_gene_transcripts=False, tpm=100.0)
    starved = I.analyze(cfg, sequences=seqs, background_gene_transcripts=False, tpm=1e-4)

    assert not rich["contrast"]["beyond_linear"], rich["contrast"]
    assert starved["contrast"]["beyond_linear"], starved["contrast"]
    assert (starved["contrast"]["min_resolvable_log2fc"]
            > rich["contrast"]["min_resolvable_log2fc"])
    for e in list(rich["groups"].values()) + [rich["contrast"]]:
        assert math.isfinite(e["min_resolvable_log2fc"])


def test_coherence_separates_a_real_group_from_a_shared_acceptor():
    seqs = _rng_seqs({"A0": ("A", 400), "A1": ("A", 380),
                      "B0": ("B", 400), "B1": ("B", 380)}, seed=9)
    groups = {"A": ["A0", "A1"], "B": ["B0", "B1"]}
    _, _, _, tracks = _system(seqs, groups)
    coherent = I.class_coherence(tracks, groups["A"])
    assert coherent["min_jaccard"] > 0.3, coherent

    seqs["X"] = "".join(random.Random(1).choice("ACGT") for _ in range(900))
    tracks["X"] = I.kmer_track(seqs["X"], W, True)
    incoherent = I.class_coherence(tracks, ["A0", "X"])
    assert incoherent["min_jaccard"] < 0.01, incoherent
    assert I.class_coherence(tracks, ["A0"])["median_jaccard"] is None


@pytest.mark.parametrize("min_log2fc,expect", [(None, "conditioning"), (0.001, "log2FC")])
def test_min_log2fc_replaces_tau_in_the_reasons(min_log2fc, expect):
    seqs = _rng_seqs({"A0": ("A", 400), "A1": ("A", 300),
                      "B0": ("B", 400), "B1": ("B", 300)}, seed=13)
    cfg = {"groups": {"A": ["A0", "A1"], "B": ["B0", "B1"]},
           "primary_comparison": ["A", "B"]}
    res = I.analyze(cfg, sequences=seqs, background_gene_transcripts=False,
                    conditioning_tau=0.0, min_log2fc=min_log2fc)
    joined = " ".join(res["reasons"])
    assert expect in joined, joined
    assert res["design"]["min_log2fc"] == min_log2fc
