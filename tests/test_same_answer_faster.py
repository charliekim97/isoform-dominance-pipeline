"""The compatibility system is built faster, and is the same system.

A ``--background-fasta`` record that shares a window is a column now, and a repeat can make
thousands of records share one: on a synthetic case, 2,364 competitors took 164 s and
3.3 GB, 55 s of it in ``compatibility_matrix``, which built a new frozenset at every
position and compared it element by element at each dictionary lookup, and most of the
rest in four singular value decompositions of one matrix.  Each function below is compared
with its 2.4.1 form, kept here as the reference, and must give the same answer exactly.
"""
import random

import numpy as np

from isoform_dominance import identifiability as I


def _kmer_track_241(seq, k, canonical=True):
    seq = seq.upper()
    if k <= 0 or len(seq) < k:
        return []
    if canonical:
        return [I.canonical_kmer(seq[i:i + k]) for i in range(len(seq) - k + 1)]
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def _compatibility_matrix_241(tracks, transcript_ids):
    sig = {}
    for tid in transcript_ids:
        for w in tracks.get(tid, []):
            sig.setdefault(w, set()).add(tid)
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


def test_kmer_track_is_the_window_by_window_fold():
    r = random.Random(11)
    for _ in range(500):
        alphabet = r.choice(["ACGT", "ACGTN", "acgtACGT", "ACGTRYN-", "AC", "ACGTÄß"])
        seq = "".join(r.choice(alphabet) for _ in range(r.randint(0, 120)))
        k = r.randint(-1, 40)
        for canonical in (True, False):
            assert I.kmer_track(seq, k, canonical) == _kmer_track_241(seq, k, canonical)
            assert I.kmers(seq, k, canonical) == set(_kmer_track_241(seq, k, canonical))


def _random_tracks(r):
    parts = ["".join(r.choice("ACGT") for _ in range(r.randint(3, 40))) for _ in range(6)]
    parts.append("A" * r.randint(5, 30))                     # a run: repeated windows
    k = r.randint(2, 7)
    tracks, ids = {}, []
    for j in range(r.randint(1, 9)):
        t = "T%d" % j
        seq = "".join(r.choice(parts) for _ in range(r.randint(0, 5)))
        tracks[t] = I.kmer_track(seq, k, r.random() < 0.5)   # may be empty: shorter than k
        ids.append(t)
    if r.random() < 0.2:
        ids.append("missing")                                # an id with no track at all
    if r.random() < 0.2:
        ids.append(ids[0])                                   # an id listed twice
    return tracks, ids


def test_compatibility_matrix_is_the_same_matrix():
    r = random.Random(12)
    shared_classes = 0
    for _ in range(800):
        tracks, ids = _random_tracks(r)
        A, classes = I.compatibility_matrix(tracks, ids)
        A0, classes0 = _compatibility_matrix_241(tracks, ids)
        assert classes == classes0
        assert A.shape == A0.shape and np.array_equal(A, A0)
        shared_classes += any(len(c) > 1 for c in classes)
    assert shared_classes > 400                              # the cases share windows


def test_one_decomposition_gives_each_functional_its_own_answer():
    r = random.Random(13)
    for _ in range(200):
        tracks, ids = _random_tracks(r)
        A, _ = I.compatibility_matrix(tracks, ids)
        if not A.size:
            continue
        svd = np.linalg.svd(A.T, full_matrices=False)[:2]
        for c in (np.ones(A.shape[1]), np.eye(A.shape[1])[0],
                  np.array([r.choice([-1.0, 0.0, 1.0]) for _ in range(A.shape[1])])):
            assert I.estimability(A, c, svd=svd) == I.estimability(A, c)
