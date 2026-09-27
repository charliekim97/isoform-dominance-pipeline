"""Gene symbol -> proposed isoform groups via the Ensembl REST API.

The hard part of isoform analysis is deciding which transcripts form a functional
group. This module clusters a gene's protein-coding transcripts by their 3' terminal-exon
splice acceptor — the alternative last exon that distinguishes functional isoform classes
(e.g. a long signalling form vs a short truncated form) — and proposes a two-group
comparison (canonical-isoform cluster vs the largest alternative cluster) that the user
reviews and renames before use.
"""
import json

from . import ensembl
from .ensembl import DEFAULT_RETRIES, DEFAULT_RETRY_WAIT

ENSEMBL = ensembl.SERVER


def _get(path, timeout=ensembl.DEFAULT_TIMEOUT, retries=DEFAULT_RETRIES,
         retry_wait=DEFAULT_RETRY_WAIT, server=None):
    return ensembl.get_json(path, server=server, timeout=timeout, retries=retries,
                            retry_wait=retry_wait)


def fetch_transcripts(gene, species="homo_sapiens", retries=DEFAULT_RETRIES,
                      retry_wait=DEFAULT_RETRY_WAIT, server=None):
    """Return {gene, species, strand, transcripts:[{id, protein_aa, terminal_acceptor, is_canonical}]}.

    The lookup is retried as described in :mod:`isoform_dominance.ensembl`.  ``server``
    is a base URL from :func:`isoform_dominance.ensembl.resolve_server`; None means the
    current release.
    """
    g = _get("/lookup/symbol/%s/%s?expand=1" % (species, gene),
             retries=retries, retry_wait=retry_wait, server=server)
    strand = g["strand"]
    canonical = (g.get("canonical_transcript") or "").split(".")[0]
    out = []
    for t in g.get("Transcript", []):
        if t.get("biotype") != "protein_coding":
            continue
        tl = t.get("Translation") or {}
        plen = tl.get("length")
        if not plen:
            continue
        exons = t["Exon"]
        if strand == 1:
            term = max(exons, key=lambda e: e["end"]); acc = term["start"]
        else:
            term = min(exons, key=lambda e: e["start"]); acc = term["end"]
        tid = t["id"].split(".")[0]
        out.append({"id": tid, "protein_aa": plen, "terminal_acceptor": acc,
                    "is_canonical": (t.get("is_canonical", 0) == 1) or (tid == canonical)})
    if not out:
        raise ValueError("No protein-coding transcripts with a translation found for %s" % gene)
    return {"gene": gene, "species": species, "strand": strand, "transcripts": out}


def cluster_by_terminal_exon(info):
    """Group transcripts by 3' terminal-exon acceptor coordinate."""
    clusters = {}
    for t in info["transcripts"]:
        clusters.setdefault(t["terminal_acceptor"], []).append(t)
    out = []
    for acc, txs in clusters.items():
        lens = sorted(x["protein_aa"] for x in txs)
        out.append({"acceptor": acc, "rep_aa": lens[len(lens) // 2], "n": len(txs),
                    "canonical": any(x["is_canonical"] for x in txs),
                    "ids": sorted(x["id"] for x in txs)})
    return sorted(out, key=lambda c: -c["n"])


#: How :func:`propose_groups` picks the alternative class; recorded in every config.
ALTERNATIVE_RULE = ("the non-canonical cluster with the most transcripts; a tie goes to the "
                    "longer representative protein, and a tie on both to the lower terminal-"
                    "acceptor coordinate")


def alternative_ties(clusters, groups, primary):
    """The clusters the alternative class was chosen over on protein length alone.

    :func:`propose_groups` takes the non-canonical cluster with the most transcripts.  When
    another cluster has as many, the choice fell to representative protein length -- a
    tie-break, not a biological criterion -- and a later release that adds one transcript
    to either cluster changes the proposed comparison.  Across six Ensembl releases (110,
    112-116), 26 to 39 of the 84 to 100 two-class proposals in a 109-gene survey were ties
    of this kind.  Of the 84 genes proposed a pair at every one of them, 30 were proposed a
    different pair at some release; 12 of the 30 had been a tie at the release before the
    first change, and in 8 -- LEPR and NTRK3 among them -- the new alternative was a
    cluster the old one had been tied with.  In all 30 it was the alternative class that
    changed, never the canonical one.
    """
    if len(primary) < 2:
        return []
    alt_ids, canon_ids = set(groups[primary[0]]), set(groups[primary[1]])
    alt = next(c for c in clusters if set(c["ids"]) == alt_ids)
    return [c for c in clusters
            if set(c["ids"]) not in (alt_ids, canon_ids) and c["n"] == alt["n"]]


def propose_groups(info):
    """The canonical cluster, the alternative (:data:`ALTERNATIVE_RULE`) and all clusters.

    Every tie is broken by content, down to the acceptor coordinate, so the proposal is a
    function of the annotation alone.  Before 2.4 a tie on both transcript count and
    protein length fell to the order in which the server listed the transcripts, which
    nothing guarantees: 249 of the 14,054 two-class genes of GENCODE 50 (1.8%) are such
    ties, and for FOXO1 and STK11 the REST order and the GTF order pick differently.
    """
    clusters = cluster_by_terminal_exon(info)
    canon = (next((c for c in clusters if c["canonical"]), None)
             or max(clusters, key=lambda c: (c["rep_aa"], -c["acceptor"])))
    others = sorted([c for c in clusters if c is not canon],
                    key=lambda c: (-c["n"], -c["rep_aa"], c["acceptor"]))
    alt = others[0] if others else None
    lbl_canon = "iso_%daa" % canon["rep_aa"]
    groups = {lbl_canon: canon["ids"]}
    primary = [lbl_canon]
    if alt:
        lbl_alt = "iso_%daa" % alt["rep_aa"]
        if lbl_alt == lbl_canon:
            lbl_alt += "_alt"
        groups[lbl_alt] = alt["ids"]
        primary = [lbl_alt, lbl_canon]  # alternative (often shorter) first
    return groups, primary, clusters


def build_config(gene, species="homo_sapiens", release=None, **retry):
    """Build a reviewable config.json dict for `gene` from Ensembl annotation.

    ``release`` is the Ensembl release to propose the groups from; None means the one
    ``rest.ensembl.org`` currently serves, and an earlier one is read from Ensembl's
    REST archive (see :func:`isoform_dominance.ensembl.resolve_server`).
    ``ensembl_release`` records the release the groups were proposed from, as the server
    reported it.  The identifiability verdict is a function of that release, so a
    config without it cannot be re-run to the same answer.
    """
    server = ensembl.resolve_server(release, **retry)
    info = fetch_transcripts(gene, species, server=server, **retry)
    release = ensembl.release_number(_get("/info/data", server=server, **retry))
    groups, primary, clusters = propose_groups(info)
    ties = alternative_ties(clusters, groups, primary)
    return {
        "gene": gene, "species": species, "ensembl_release": release,
        "groups": groups, "primary_comparison": primary,
        "_proposal": {"alternative_rule": ALTERNATIVE_RULE,
                      "tied_with": [{"terminal_acceptor": c["acceptor"],
                                     "rep_protein_aa": c["rep_aa"],
                                     "n_transcripts": c["n"],
                                     "transcripts": c["ids"]} for c in ties]},
        "_proposed": ("Auto-proposed by `isoform-dominance annotate`. Groups = protein-coding "
                      "transcripts sharing a 3' terminal-exon splice acceptor (isoform-defining "
                      "alternative last exon). REVIEW and rename to functional names "
                      "(e.g. short/long) before use; smaller clusters are listed under _clusters."),
        "_clusters": [{"terminal_acceptor": c["acceptor"], "rep_protein_aa": c["rep_aa"],
                       "n_transcripts": c["n"], "contains_canonical": c["canonical"],
                       "transcripts": c["ids"]} for c in clusters],
    }


def run(gene, out, species="homo_sapiens", release=None, **retry):
    cfg = build_config(gene, species, release=release, **retry)
    with open(out, "w") as f:
        json.dump(cfg, f, indent=2)
    return cfg
