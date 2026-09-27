"""Shared IO helpers: config and sample-map loading."""
import csv
import json


class InputError(ValueError):
    """A file named on the command line that cannot be used as it is."""


def load_json(path, what):
    """Parse the JSON file ``path``; :class:`InputError` naming ``what`` if it is not JSON."""
    with open(path) as f:
        try:
            return json.load(f)
        except json.JSONDecodeError as e:
            raise InputError("%s %s is not valid JSON (%s)" % (what, path, e)) from e


def load_config(path, need_groups=True):
    """A config written by ``annotate`` or by hand: a JSON object, with a ``groups`` object
    mapping group names to transcript ids unless ``need_groups`` is false (``qc`` reads
    only ``contamination_qc``)."""
    cfg = load_json(path, "config")
    if not isinstance(cfg, dict):
        raise InputError("config %s is not a JSON object" % path)
    if need_groups and not isinstance(cfg.get("groups"), dict):
        raise InputError("config %s has no \"groups\" object mapping group names to "
                         "transcript ids; `annotate` writes one" % path)
    return cfg



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

