"""Structural invariants of the compatibility system.

These are theorems about the matrix, not observations about any gene, so they are
stated as tests: if a change to ``compatibility_matrix`` or ``estimability`` breaks
one, the change is wrong.

The theorem the first three tests rest on: every column of ``A`` sums to 1, so
``1_rows' A = 1_cols'``; if ``A v = 0`` then ``1_cols' v = 0``.  Every null direction
therefore sums to zero over *all* transcripts, and the gene total is always estimable.
Its precondition -- every transcript at least ``window`` long -- is the subject of the
last test, which is what happens when the precondition fails.
"""
import random

import numpy as np
import pytest

from isoform_dominance import identifiability as I

W = 31


def _system(seqs, ids=None, window=W):
    ids = sorted(seqs) if ids is None else ids
    tracks = {t: I.kmer_track(seqs[t], window, True) for t in ids}
    out = I.compatibility_matrix(tracks, ids)
    return (out[0] if isinstance(out, tuple) else out), ids


def _null_space(A):
    _, s, vt = np.linalg.svd(A, full_matrices=True)
    tol = max(A.shape) * (s[0] if s.size else 0.0) * 1e-10
    return vt[int((s > tol).sum()):]


def _seqs(spec, seed=0):
    rng = random.Random(seed)

    def draw(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    core = draw(500)
    out = {}
    for name, kind in spec.items():
        out[name] = draw(400) if kind == "random" else core + draw(120)
    return out, core


@pytest.mark.parametrize("spec,extra,min_nullity", [
    ({"T%d" % i: "random" for i in range(8)}, None, 0),
    ({"T%d" % i: "shared" for i in range(8)}, None, 0),
    ({"T%d" % i: "shared" for i in range(6)}, "duplicate", 1),
    ({"T%d" % i: "shared" for i in range(4)}, "identical", 3),
])
def test_the_gene_total_is_always_estimable(spec, extra, min_nullity):
    """Forced by column-stochasticity; a failure means the columns stopped summing to 1."""
    seqs, core = _seqs(spec)
    if extra == "duplicate":
        seqs["DUP"] = seqs["T0"]
    elif extra == "identical":
        seqs = dict.fromkeys(seqs, core[:200])
    A, ids = _system(seqs)
    assert np.allclose(A.sum(axis=0), 1.0), "columns must sum to 1"
    null = _null_space(A)
    assert null.shape[0] >= min_nullity, (
        "this case is meant to exercise a non-trivial null space; with none, the "
        "per-direction assertion below is vacuous")
    for v in null:
        assert abs(v.sum()) < 1e-9, "a null direction did not sum to zero over all transcripts"
    e = I.estimability(A, np.ones(len(ids)))
    assert e["estimable"] and e["residual"] < 1e-9


def test_a_class_total_is_estimable_exactly_when_no_null_direction_crosses_it():
    """Within-class null directions are harmless; one that does not cancel is fatal.

    Two identical transcripts in the same class leave the class total exactly
    recoverable even though neither abundance is.  Move one of them into the other
    class and the same null direction now has weight 1 on one side and 0 on the other,
    so the class indicator leaves the row space.
    """
    seqs, _ = _seqs({"A0": "shared", "A1": "shared", "B0": "random", "B1": "random"})
    seqs["A2"] = seqs["A0"]                      # exact duplicate of A0

    for cls_of_a2, expected in (("A", True), ("B", False)):
        groups = {"A": ["A0", "A1"], "B": ["B0", "B1"]}
        groups[cls_of_a2].append("A2")
        ids = groups["A"] + groups["B"]
        A, ids = _system(seqs, ids)
        c = np.array([1.0 if t in groups["A"] else 0.0 for t in ids])

        crossing = [v for v in _null_space(A) if abs(v @ c) > 1e-9]
        e = I.estimability(A, c)
        assert e["estimable"] is expected
        assert bool(crossing) is (not expected), (
            "estimability must agree with whether a null direction crosses the class")
        assert I.estimability(A, np.ones(len(ids)))["estimable"], "gene total, either way"


def test_a_transcript_shorter_than_the_window_breaks_column_stochasticity():
    """The precondition of the theorem above, and what its failure costs.

    Such a transcript has no windows, so its column is all zero, ``A`` is no longer
    column-stochastic, and the guarantee is void -- the *gene* total stops being
    estimable.  Callers should drop or reject these transcripts rather than report the
    resulting ``not_identifiable``.
    """
    seqs, _ = _seqs({"T%d" % i: "shared" for i in range(4)})
    seqs["SHORT"] = "ACGT" * 3                                    # 12 nt < window
    assert I.kmer_track(seqs["SHORT"], W, True) == []

    A, ids = _system(seqs)
    col = A[:, ids.index("SHORT")]
    assert not col.any(), "a transcript with no windows must give an all-zero column"
    assert not np.allclose(A.sum(axis=0), 1.0)

    e = I.estimability(A, np.ones(len(ids)))
    assert not e["estimable"], (
        "documents the cost of the broken precondition; if a guard is added upstream "
        "so that this can no longer happen, delete this test with it")
