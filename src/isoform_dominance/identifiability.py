"""Are two isoform groups distinguishable by short reads?

Short-read quantifiers (e.g. Salmon) can only apportion isoforms using sequence that is
UNIQUE to one isoform group. If a group shares all of its sequence with another group, no
short read can be assigned to it unambiguously and its abundance is not identifiable.

This module fetches transcript cDNA (Ensembl, or supplied sequences), builds per-group
k-mer sets, and reports each group's group-unique k-mer count. A group with zero unique
k-mers is flagged as NOT distinguishable by short reads — an honest guardrail that most
isoform analyses skip.
"""
import json
import urllib.request

ENSEMBL = "https://rest.ensembl.org"


def _get_text(path, timeout=30):
    req = urllib.request.Request(ENSEMBL + path, headers={"Content-Type": "text/plain"})
    return urllib.request.urlopen(req, timeout=timeout).read().decode().strip()


def fetch_cdna(transcript_id):
    return _get_text("/sequence/id/%s?type=cdna" % transcript_id.split(".")[0])


def kmers(seq, k):
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def analyze(config, k=31, sequences=None):
    """sequences: optional {transcript_id: cdna}. If None, fetched from Ensembl.
    Returns {group: {n_unique_kmers, n_transcripts, distinguishable}} + summary."""
    groups = config["groups"]
    needed = {t.split(".")[0] for ids in groups.values() for t in ids}
    seqs = dict(sequences or {})
    for tid in needed:
        if tid not in seqs:
            seqs[tid] = fetch_cdna(tid)

    group_kmers = {}
    for g, ids in groups.items():
        ks = set()
        for t in ids:
            ks |= kmers(seqs[t.split(".")[0]], k)
        group_kmers[g] = ks

    report = {}
    for g, ks in group_kmers.items():
        others = set()
        for g2, ks2 in group_kmers.items():
            if g2 != g:
                others |= ks2
        uniq = ks - others
        report[g] = {"n_unique_kmers": len(uniq), "n_transcripts": len(groups[g]),
                     "distinguishable": len(uniq) > 0}

    pc = config.get("primary_comparison", list(groups)[:2])
    primary_ok = all(report[g]["distinguishable"] for g in pc if g in report)
    return {"k": k, "groups": report, "primary_comparison": pc,
            "primary_distinguishable": primary_ok}


def run(config, k=31, sequences=None):
    res = analyze(config, k=k, sequences=sequences)
    return res
