"""`annotate` against a recorded Ensembl payload.

The Ensembl-parsing path is the one part of the package that cannot be exercised by
supplying inputs directly, so it was previously untested and rotted silently whenever
the REST schema shifted.  The fixture below is a trimmed but structurally faithful
``/lookup/symbol/...?expand=1`` response modelled on human *LEPR*: a minus-strand gene
whose canonical long isoform (LepRb) uses a distal terminal exon and whose short
isoforms (LepRa) share an earlier one, plus the non-coding and untranslated transcripts
that a real payload carries and that the clustering has to discard.
"""
import json

import pytest

from isoform_dominance import annotate

# minus-strand gene: the terminal exon is the one with the smallest start, and the
# acceptor of interest is its `end`
LEPR_LIKE = {
    "id": "ENSG00000116678",
    "display_name": "LEPR",
    "strand": -1,
    "canonical_transcript": "ENST00000349533.11",
    "Transcript": [
        # --- long class: canonical, distal terminal exon (acceptor 65_580_000) ----
        {"id": "ENST00000349533", "biotype": "protein_coding", "is_canonical": 1,
         "Translation": {"length": 1165},
         "Exon": [{"start": 65_580_000, "end": 65_583_000},
                  {"start": 65_600_000, "end": 65_600_400},
                  {"start": 65_650_000, "end": 65_650_300}]},
        {"id": "ENST00000344610", "biotype": "protein_coding", "is_canonical": 0,
         "Translation": {"length": 1162},
         "Exon": [{"start": 65_579_500, "end": 65_583_000},
                  {"start": 65_600_000, "end": 65_600_400}]},
        # --- short class: proximal terminal exon (acceptor 65_571_200) ------------
        # transcripts of one class share the terminal-exon splice acceptor and differ
        # only in where the 3' UTR ends -- on the minus strand that is the exon `start`
        {"id": "ENST00000371060", "biotype": "protein_coding", "is_canonical": 0,
         "Translation": {"length": 896},
         "Exon": [{"start": 65_570_000, "end": 65_571_200},
                  {"start": 65_600_000, "end": 65_600_400}]},
        {"id": "ENST00000616738", "biotype": "protein_coding", "is_canonical": 0,
         "Translation": {"length": 896},
         "Exon": [{"start": 65_569_400, "end": 65_571_200},
                  {"start": 65_600_000, "end": 65_600_400}]},
        {"id": "ENST00000371058", "biotype": "protein_coding", "is_canonical": 0,
         "Translation": {"length": 892},
         "Exon": [{"start": 65_568_800, "end": 65_571_200}]},
        # --- must be discarded ----------------------------------------------------
        {"id": "ENST00000900001", "biotype": "retained_intron",
         "Exon": [{"start": 65_560_000, "end": 65_562_000}]},
        {"id": "ENST00000900002", "biotype": "protein_coding", "Translation": None,
         "Exon": [{"start": 65_555_000, "end": 65_556_000}]},
        {"id": "ENST00000900003", "biotype": "processed_transcript",
         "Exon": [{"start": 65_550_000, "end": 65_551_000}]},
    ],
}


@pytest.fixture
def offline_ensembl(monkeypatch):
    calls = []

    def fake_get(path, timeout=30):
        calls.append(path)
        return LEPR_LIKE

    monkeypatch.setattr(annotate, "_get", fake_get)
    return calls


# --------------------------------------------------------------------------- #
def test_fetch_transcripts_keeps_only_translated_protein_coding(offline_ensembl):
    info = annotate.fetch_transcripts("LEPR")
    ids = {t["id"] for t in info["transcripts"]}
    assert ids == {"ENST00000349533", "ENST00000344610", "ENST00000371060",
                   "ENST00000616738", "ENST00000371058"}
    assert info["strand"] == -1
    assert offline_ensembl == ["/lookup/symbol/homo_sapiens/LEPR?expand=1"]


def test_fetch_transcripts_uses_the_minus_strand_terminal_exon(offline_ensembl):
    info = annotate.fetch_transcripts("LEPR")
    by_id = {t["id"]: t for t in info["transcripts"]}
    # on the minus strand the terminal exon is the leftmost, and its acceptor is `end`
    assert by_id["ENST00000349533"]["terminal_acceptor"] == 65_583_000
    assert by_id["ENST00000371060"]["terminal_acceptor"] == 65_571_200
    assert by_id["ENST00000616738"]["terminal_acceptor"] == 65_571_200


def test_plus_strand_uses_the_rightmost_exon_start(monkeypatch):
    payload = json.loads(json.dumps(LEPR_LIKE))
    payload["strand"] = 1
    monkeypatch.setattr(annotate, "_get", lambda path, timeout=30: payload)
    info = annotate.fetch_transcripts("LEPR")
    by_id = {t["id"]: t for t in info["transcripts"]}
    assert by_id["ENST00000349533"]["terminal_acceptor"] == 65_650_000


def test_canonical_flag_survives_the_version_suffix(offline_ensembl):
    # canonical_transcript is "ENST00000349533.11"; the transcript id has no version
    info = annotate.fetch_transcripts("LEPR")
    canon = [t["id"] for t in info["transcripts"] if t["is_canonical"]]
    assert canon == ["ENST00000349533"]


def test_clusters_are_ordered_by_size(offline_ensembl):
    info = annotate.fetch_transcripts("LEPR")
    clusters = annotate.cluster_by_terminal_exon(info)
    assert [c["n"] for c in clusters] == sorted([c["n"] for c in clusters], reverse=True)
    assert sum(c["n"] for c in clusters) == 5
    assert any(c["canonical"] for c in clusters)


def test_propose_groups_pairs_the_canonical_cluster_with_the_largest_alternative(offline_ensembl):
    info = annotate.fetch_transcripts("LEPR")
    groups, primary, clusters = annotate.propose_groups(info)
    assert len(groups) == 2
    canon_label = next(lab for lab, ids in groups.items()
                       if "ENST00000349533" in ids)
    alt_label = next(lab for lab in groups if lab != canon_label)
    # the alternative (shorter) class is named first, so it is the numerator
    assert primary == [alt_label, canon_label]
    assert set(groups[alt_label]) == {"ENST00000371060", "ENST00000616738",
                                      "ENST00000371058"}
    assert canon_label == "iso_1165aa"


def test_build_config_is_reviewable_and_records_every_cluster(offline_ensembl):
    cfg = annotate.build_config("LEPR")
    assert cfg["gene"] == "LEPR"
    assert "REVIEW" in cfg["_proposed"]
    assert len(cfg["_clusters"]) >= 2
    accs = {c["terminal_acceptor"] for c in cfg["_clusters"]}
    assert {65_583_000, 65_571_200} <= accs
    # every clustered transcript is accounted for
    assert sum(c["n_transcripts"] for c in cfg["_clusters"]) == 5


def test_run_writes_the_config(offline_ensembl, tmp_path):
    out = tmp_path / "cfg.json"
    cfg = annotate.run("LEPR", str(out))
    assert json.load(open(out)) == cfg
    assert cfg["primary_comparison"][0] in cfg["groups"]


def test_no_translated_transcripts_raises(monkeypatch):
    payload = {"strand": 1, "canonical_transcript": "",
               "Transcript": [{"id": "T1", "biotype": "lncRNA",
                               "Exon": [{"start": 1, "end": 100}]}]}
    monkeypatch.setattr(annotate, "_get", lambda path, timeout=30: payload)
    with pytest.raises(ValueError, match="No protein-coding transcripts"):
        annotate.fetch_transcripts("NOPE")


def test_single_cluster_gene_proposes_one_group(monkeypatch):
    """A gene with no alternative terminal exon must not fabricate a comparison."""
    payload = {"strand": 1, "canonical_transcript": "T1",
               "Transcript": [
                   {"id": "T1", "biotype": "protein_coding", "is_canonical": 1,
                    "Translation": {"length": 300},
                    "Exon": [{"start": 1, "end": 100}, {"start": 200, "end": 300}]},
                   {"id": "T2", "biotype": "protein_coding", "is_canonical": 0,
                    "Translation": {"length": 290},
                    "Exon": [{"start": 1, "end": 100}, {"start": 200, "end": 290}]}]}
    monkeypatch.setattr(annotate, "_get", lambda path, timeout=30: payload)
    groups, primary, _ = annotate.propose_groups(annotate.fetch_transcripts("X"))
    assert len(groups) == 1
    assert len(primary) == 1


def test_duplicate_labels_are_disambiguated(monkeypatch):
    """Two clusters whose median protein length matches must not collide on one label."""
    payload = {"strand": 1, "canonical_transcript": "T1",
               "Transcript": [
                   {"id": "T1", "biotype": "protein_coding", "is_canonical": 1,
                    "Translation": {"length": 500},
                    "Exon": [{"start": 1, "end": 100}, {"start": 400, "end": 500}]},
                   {"id": "T2", "biotype": "protein_coding", "is_canonical": 0,
                    "Translation": {"length": 500},
                    "Exon": [{"start": 1, "end": 100}, {"start": 700, "end": 800}]}]}
    monkeypatch.setattr(annotate, "_get", lambda path, timeout=30: payload)
    groups, primary, _ = annotate.propose_groups(annotate.fetch_transcripts("X"))
    assert len(groups) == 2                       # not one label overwriting the other
    assert len(set(groups)) == 2


def test_clustering_is_brittle_to_a_shifted_acceptor(monkeypatch):
    """Known limitation, pinned so a fix is a visible behaviour change.

    Clusters are keyed on the exact terminal-exon acceptor coordinate, so two
    transcripts of the same functional class whose annotated acceptor differs by even
    one base land in different clusters.  Annotation sets do contain such shifts, and
    the effect is a silently mis-specified group rather than an error.  An acceptor
    tolerance is the intended fix; until then this test documents the behaviour.
    """
    payload = json.loads(json.dumps(LEPR_LIKE))
    for t in payload["Transcript"]:
        if t["id"] == "ENST00000616738":
            t["Exon"][0]["end"] = 65_571_199          # one base off
    monkeypatch.setattr(annotate, "_get", lambda path, timeout=30: payload)
    clusters = annotate.cluster_by_terminal_exon(annotate.fetch_transcripts("LEPR"))
    accs = sorted(c["acceptor"] for c in clusters)
    assert 65_571_199 in accs and 65_571_200 in accs   # split, not merged
