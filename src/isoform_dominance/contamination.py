"""Contamination control: does a target isoform's signal track a contaminating cell type?

Matplotlib is imported lazily inside ``run`` to keep CLI startup fast for
subcommands that do not plot.
"""
import csv, os, math
import numpy as np
from scipy.stats import spearmanr

from .io import InputError

PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]


def _log2p1(x):
    return math.log2(x + 1.0)


def load_markers(path, tissue, contaminant):
    out = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        cols = set(reader.fieldnames or [])
        have_t = [g for g in tissue if g in cols]
        have_c = [g for g in contaminant if g in cols]
        if not have_t or not have_c:
            raise InputError(
                "%s is missing marker columns: present tissue=%s, contaminant=%s "
                "(need at least one column from each panel)." % (path, have_t, have_c))
        for r in reader:
            t = np.mean([_log2p1(float(r[g])) for g in have_t])
            c = np.mean([_log2p1(float(r[g])) for g in have_c])
            out[r["donor"]] = (t, c, (c / t if t > 0 else float("nan")))
    return out


def load_target(path, target_group):
    col = "%s_TPM" % target_group
    with open(path) as f:
        reader = csv.DictReader(f)
        if col not in (reader.fieldnames or []):
            raise InputError("%s has no column %r (target_group=%s)." % (path, col, target_group))
        return {r["donor"]: float(r[col]) for r in reader}


_QC_SHAPE = ("add it with a 'target_group' and 'marker_panels' ({'tissue': [...], "
             "'contaminant': [...]}) to run qc; `annotate` does not write one")


def qc_section(config):
    """``(target_group, tissue markers, contaminant markers)`` from the config's
    ``contamination_qc``; :class:`InputError` saying what is missing."""
    qc = config.get("contamination_qc")
    if not isinstance(qc, dict):
        raise InputError("config has no 'contamination_qc' section; " + _QC_SHAPE)
    if not isinstance(qc.get("target_group"), str):
        raise InputError("config's contamination_qc has no 'target_group'; " + _QC_SHAPE)
    panels = qc.get("marker_panels")
    if not isinstance(panels, dict):
        raise InputError("config's contamination_qc has no 'marker_panels'; " + _QC_SHAPE)
    for key in ("tissue", "contaminant"):
        if not (isinstance(panels.get(key), list) and panels[key]):
            raise InputError("config's contamination_qc.marker_panels has no '%s' list of "
                             "marker genes; %s" % (key, _QC_SHAPE))
    return qc["target_group"], panels["tissue"], panels["contaminant"]


def run(config, markers, targets, out):
    """markers/targets: {cohort: path}. Writes <out>.{png,pdf,svg}+_scores.csv. Returns rows."""
    import matplotlib as mpl; mpl.use("Agg")
    import matplotlib.pyplot as plt
    tg, tissue, contam = qc_section(config)
    names = list(markers)
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 9, "pdf.fonttype": 42, "svg.fonttype": "none"})
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig, axes = plt.subplots(1, len(names), figsize=(3.6 * len(names), 3.4), squeeze=False)
    srows = []
    for i, name in enumerate(names):
        m = load_markers(markers[name], tissue, contam)
        t = load_target(targets[name], tg)
        donors = [d for d in m if d in t]
        if len(donors) < 3:
            raise InputError(
                "cohort %s: only %d donor(s) overlap between markers and target tables; "
                "need >= 3 for a Spearman correlation." % (name, len(donors)))
        cs = np.array([m[d][1] for d in donors])
        ratio = np.array([m[d][2] for d in donors])
        tv = np.array([t[d] for d in donors])
        rho, p = spearmanr(cs, tv)
        srows.append((name, len(donors), rho, p, float(np.nanmedian(ratio))))
        ax = axes[0, i]; c = PALETTE[i % len(PALETTE)]
        ax.scatter(cs, tv, s=36, c=c, edgecolor="white", lw=0.5, zorder=3)
        ax.set_xlabel("contamination score\n[mean log2(TPM+1), contaminant]")
        ax.set_ylabel("%s (%s) TPM" % (config.get("gene", "target"), tg))
        ax.set_title("%s (n=%d)" % (name, len(donors)), fontsize=9.5)
        ax.text(0.03, 0.97, "rho=%.2f\nP=%.2f" % (rho, p), transform=ax.transAxes, va="top",
                fontsize=8, bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#ccc", "lw": 0.6})
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("%s (%s) vs contamination - no positive dependence = genuine signal"
                 % (config.get("gene", "target"), tg), fontsize=10, fontweight="bold", y=1.04)
    for ext in ("png", "pdf", "svg"):
        fig.savefig("%s.%s" % (out, ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    with open(out + "_scores.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["cohort", "n", "spearman_rho", "spearman_P", "median_contam_tissue_ratio"])
        for r in srows:
            w.writerow([r[0], r[1], "%.3f" % r[2], "%.4g" % r[3], "%.4f" % r[4]])
    return srows
