"""Gene symbol -> proposed isoform groups via the Ensembl REST API.

The hard part of isoform analysis is deciding which transcripts form a functional
group. This module clusters a gene's protein-coding transcripts by their 3' terminal-exon
splice acceptor — the alternative last exon that distinguishes functional isoform classes
(e.g. a long signalling form vs a short truncated form) — and proposes a two-group
comparison (canonical-isoform cluster vs the largest alternative cluster) that the user
reviews and renames before use.
"""
import json
import urllib.request

ENSEMBL = "https://rest.ensembl.org"


def _get(path, timeout=30):
    req = urllib.request.Request(ENSEMBL + path, headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=timeout))


def fetch_transcripts(gene, species="homo_sapiens"):
    """Return {gene, species, strand, transcripts:[{id, protein_aa, terminal_acceptor, is_canonical}]}."""
    g = _get("/lookup/symbol/%s/%s?expand=1" % (species, gene))
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


def propose_groups(info):
    clusters = cluster_by_terminal_exon(info)
    canon = next((c for c in clusters if c["canonical"]), None) or max(clusters, key=lambda c: c["rep_aa"])
    others = sorted([c for c in clusters if c is not canon], key=lambda c: (-c["n"], -c["rep_aa"]))
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


def build_config(gene, species="homo_sapiens"):
    """Build a reviewable config.json dict for `gene` from live Ensembl annotation."""
    info = fetch_transcripts(gene, species)
    groups, primary, clusters = propose_groups(info)
    return {
        "gene": gene, "species": species, "reference": "Ensembl REST (live annotation)",
        "groups": groups, "primary_comparison": primary,
        "_proposed": ("Auto-proposed by `isoform-dominance annotate`. Groups = protein-coding "
                      "transcripts sharing a 3' terminal-exon splice acceptor (isoform-defining "
                      "alternative last exon). REVIEW and rename to functional names "
                      "(e.g. short/long) before use; smaller clusters are listed under _clusters."),
        "_clusters": [{"terminal_acceptor": c["acceptor"], "rep_protein_aa": c["rep_aa"],
                       "n_transcripts": c["n"], "contains_canonical": c["canonical"],
                       "transcripts": c["ids"]} for c in clusters],
    }


def run(gene, out, species="homo_sapiens"):
    cfg = build_config(gene, species)
    with open(out, "w") as f:
        json.dump(cfg, f, indent=2)
    return cfg
