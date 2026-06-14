"""Contamination control: does a target isoform's signal track a contaminating cell type?"""
import csv, os, math
import numpy as np
from scipy.stats import spearmanr
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]


def _log2p1(x):
    return math.log2(x + 1.0)


def load_markers(path, tissue, contaminant):
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            t = np.mean([_log2p1(float(r[g])) for g in tissue if g in r])
            c = np.mean([_log2p1(float(r[g])) for g in contaminant if g in r])
            out[r["donor"]] = (t, c, (c / t if t > 0 else float("nan")))
    return out


def load_target(path, target_group):
    col = "%s_TPM" % target_group
    with open(path) as f:
        return {r["donor"]: float(r[col]) for r in csv.DictReader(f)}


def run(config, markers, targets, out):
    """markers/targets: {cohort: path}. Writes <out>.{png,pdf,svg}+_scores.csv. Returns rows."""
    qc = config["contamination_qc"]
    tg = qc["target_group"]
    tissue = qc["marker_panels"]["tissue"]; contam = qc["marker_panels"]["contaminant"]
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
                fontsize=8, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ccc", lw=0.6))
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
