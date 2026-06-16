"""Shared IO helpers: config and sample-map loading."""
import json, csv


def load_config(path):
    with open(path) as f:
        return json.load(f)


def transcript_to_group(groups):
    """{group: [ENST,...]} -> {ENST(no version): group}."""
    m = {}
    for g, txs in groups.items():
        for t in txs:
            m[t.split(".")[0]] = g
    return m


def load_sample_map(path):
    """CSV with columns donor,condition[,SRR] -> {donor: condition}."""
    cond = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            cond[r["donor"]] = r.get("condition", "NA")
    return cond


def primary_pair(config):
    """The two groups named in ``primary_comparison`` (for the paired comparison)."""
    groups = config["groups"]
    pc = config.get("primary_comparison", list(groups)[:2])
    if len(pc) < 2:
        raise ValueError(
            "primary_comparison must name two isoform groups for a paired "
            "comparison; got %r. `annotate` proposes a pair only when an "
            "alternative terminal-exon cluster exists; edit the config to define "
            "two groups." % (pc,))
    return pc[0], pc[1]
