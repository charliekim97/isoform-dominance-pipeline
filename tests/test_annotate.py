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


def _tied_info():
    # canonical 1000aa (acceptor 100); two alternatives of one transcript each
    return {"gene": "X", "species": "homo_sapiens", "strand": 1, "transcripts": [
        {"id": "T1", "protein_aa": 1000, "terminal_acceptor": 100, "is_canonical": True},
        {"id": "T2", "protein_aa": 900, "terminal_acceptor": 200, "is_canonical": False},
        {"id": "T3", "protein_aa": 800, "terminal_acceptor": 300, "is_canonical": False},
    ]}


def test_a_clear_alternative_has_no_ties():
    groups, primary, clusters = annotate.propose_groups(_info())
    assert annotate.alternative_ties(clusters, groups, primary) == []


def test_a_tie_on_transcript_count_is_reported_with_the_runner_up():
    groups, primary, clusters = annotate.propose_groups(_tied_info())
    assert primary == ["iso_900aa", "iso_1000aa"]           # the longer protein won
    ties = annotate.alternative_ties(clusters, groups, primary)
    assert [(c["acceptor"], c["rep_aa"], c["n"]) for c in ties] == [(300, 800, 1)]


def test_a_single_cluster_has_nothing_to_tie_with():
    info = _info()
    info["transcripts"] = info["transcripts"][:1]
    groups, primary, clusters = annotate.propose_groups(info)
    assert annotate.alternative_ties(clusters, groups, primary) == []


def _full_tie_info():
    # FOXO1 at release 116, reduced: the canonical cluster and two alternatives of one
    # 655 aa transcript each -- tied on transcript count and on protein length
    return [
        {"id": "C1", "protein_aa": 655, "terminal_acceptor": 40559034, "is_canonical": True},
        {"id": "C2", "protein_aa": 655, "terminal_acceptor": 40559034, "is_canonical": False},
        {"id": "X1", "protein_aa": 655, "terminal_acceptor": 40560860, "is_canonical": False},
        {"id": "Y1", "protein_aa": 655, "terminal_acceptor": 40558509, "is_canonical": False},
    ]


def test_a_full_tie_is_broken_by_content_not_by_the_order_transcripts_arrive_in():
    # before 2.4 the alternative among fully tied clusters was whichever the server listed
    # first; REST and the GTF list FOXO1's in different orders and picked differently
    import itertools
    picks = set()
    for order in itertools.permutations(_full_tie_info()):
        info = {"gene": "FOXO1", "species": "homo_sapiens", "strand": -1,
                "transcripts": list(order)}
        groups, primary, _ = annotate.propose_groups(info)
        picks.add(tuple(groups[primary[0]]))
    assert picks == {("Y1",)}           # the lower acceptor, 40558509, in every order


def test_a_missing_canonical_falls_back_by_content_too():
    import itertools
    rows = [{"id": "A", "protein_aa": 700, "terminal_acceptor": 500, "is_canonical": False},
            {"id": "B", "protein_aa": 700, "terminal_acceptor": 300, "is_canonical": False},
            {"id": "C", "protein_aa": 400, "terminal_acceptor": 900, "is_canonical": False}]
    canons = set()
    for order in itertools.permutations(rows):
        groups, primary, _ = annotate.propose_groups(
            {"gene": "X", "species": "homo_sapiens", "strand": 1, "transcripts": list(order)})
        canons.add(tuple(groups[primary[-1]]))
    assert canons == {("B",)}           # longest protein, then the lower acceptor


def test_the_recorded_rule_names_every_tie_break():
    assert "lower terminal-acceptor coordinate" in annotate.ALTERNATIVE_RULE
