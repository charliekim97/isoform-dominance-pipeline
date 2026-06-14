"""Shared IO helpers: config and sample-map loading."""
import json, csv, os


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
    groups = config["groups"]
    pc = config.get("primary_comparison", list(groups)[:2])
    return pc[0], pc[1]
