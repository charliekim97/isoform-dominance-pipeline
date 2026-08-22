"""Short-read distinguishability (offline; supplied sequences)."""
import pytest

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


def test_rejects_primary_comparison_naming_an_unknown_group():
    # The `annotate` workflow tells users to rename the proposed groups; renaming
    # `groups` but not `primary_comparison` used to report primary_distinguishable
    # True (the unknown label was silently skipped), passing the CLI's exit-0 gate.
    seqs = {"A1": SHARED + "GGGGGGGGCATCAT", "B1": SHARED + "TTTTTTTTAGAGAG"}
    cfg = {"groups": {"A": ["A1"], "B": ["B1"]}, "primary_comparison": ["A", "Bee"]}
    with pytest.raises(ValueError, match="not in config"):
        identifiability.analyze(cfg, k=8, sequences=seqs)


def test_rejects_single_group_config():
    # One group means no comparison at all: `others` is empty, so the lone group
    # is trivially "distinguishable" and the summary used to read True.
    seqs = {"S1": SHARED + "GGGGGGGG"}
    cfg = {"groups": {"only": ["S1"]}, "primary_comparison": ["only"]}
    with pytest.raises(ValueError, match="two isoform groups"):
        identifiability.analyze(cfg, k=8, sequences=seqs)
