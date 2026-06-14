"""Short-read distinguishability (offline; supplied sequences)."""
from isoform_dominance import identifiability

SHARED = "ACGT" * 20  # 80 bp shared backbone


def test_distinguishable_when_each_group_has_unique_region():
    seqs = {"A1": SHARED + "GGGGGGGGCATCAT", "B1": SHARED + "TTTTTTTTAGAGAG"}
    cfg = {"groups": {"A": ["A1"], "B": ["B1"]}, "primary_comparison": ["A", "B"]}
    res = identifiability.analyze(cfg, k=8, sequences=seqs)
    assert res["primary_distinguishable"] is True
    assert res["groups"]["A"]["n_unique_kmers"] > 0
    assert res["groups"]["B"]["n_unique_kmers"] > 0


def test_flags_group_with_no_unique_sequence():
    # 'sub' sequence is fully contained in 'super' -> sub has no unique k-mers
    sub = SHARED + "GGGGGGGG"
    sup = sub + "TTTTTTTTTT"
    seqs = {"S1": sub, "P1": sup}
    cfg = {"groups": {"sub": ["S1"], "super": ["P1"]}, "primary_comparison": ["sub", "super"]}
    res = identifiability.analyze(cfg, k=8, sequences=seqs)
    assert res["groups"]["sub"]["n_unique_kmers"] == 0
    assert res["groups"]["sub"]["distinguishable"] is False
    assert res["primary_distinguishable"] is False
