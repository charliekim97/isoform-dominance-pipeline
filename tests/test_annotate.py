"""Isoform-grouping logic (offline; no Ensembl call)."""
from isoform_dominance import annotate


def _info():
    # canonical 1000aa (acceptor 100); two 900aa share acceptor 200; one 500aa acceptor 300
    return {"gene": "X", "species": "homo_sapiens", "strand": 1, "transcripts": [
        {"id": "T1", "protein_aa": 1000, "terminal_acceptor": 100, "is_canonical": True},
        {"id": "T2", "protein_aa": 900, "terminal_acceptor": 200, "is_canonical": False},
        {"id": "T3", "protein_aa": 900, "terminal_acceptor": 200, "is_canonical": False},
        {"id": "T4", "protein_aa": 500, "terminal_acceptor": 300, "is_canonical": False},
    ]}


def test_clusters_by_terminal_exon():
    clusters = annotate.cluster_by_terminal_exon(_info())
    accs = {c["acceptor"]: c["n"] for c in clusters}
    assert accs == {200: 2, 100: 1, 300: 1}


def test_proposes_canonical_vs_largest_alt():
    groups, primary, _ = annotate.propose_groups(_info())
    # canonical cluster -> iso_1000aa; largest alternative (2 tx) -> iso_900aa
    assert groups["iso_1000aa"] == ["T1"]
    assert groups["iso_900aa"] == ["T2", "T3"]
    # alternative group is listed first in the comparison
    assert primary == ["iso_900aa", "iso_1000aa"]
