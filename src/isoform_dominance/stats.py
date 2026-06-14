"""Paired isoform-group statistics + figure."""
import csv, os
import numpy as np
from scipy.stats import wilcoxon
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from .io import primary_pair

PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]


def load_perdonor(path, condition, gA, gB):
    don, A, B = [], [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            if condition not in (None, "", "all") and r["condition"] != condition:
                continue
            don.append(r["donor"]); A.append(float(r["%s_TPM" % gA])); B.append(float(r["%s_TPM" % gB]))
    return don, np.array(A), np.array(B)


def paired_stat(A, B):
    n = len(A); ngt = int(np.sum(A > B))
    try:
        p = wilcoxon(A, B, alternative="two-sided").pvalue
    except Exception:
        p = float("nan")
    with np.errstate(divide="ignore", invalid="ignore"):
        fold = float(np.median(A / B))
    return n, ngt, p, fold


def run(config, condition, cohorts, out):
    """cohorts: {name: perdonor.csv}. Writes <out>.{png,pdf,svg} + <out>_stats.csv. Returns summary dict."""
    gene = config.get("gene", "gene")
    gA, gB = primary_pair(config)
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 9, "pdf.fonttype": 42, "svg.fonttype": "none"})
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    names = list(cohorts)
    fig, axes = plt.subplots(1, len(names), figsize=(3.6 * len(names), 3.6), squeeze=False)
    allA, allB, statrows = [], [], []
    for i, name in enumerate(names):
        don, A, B = load_perdonor(cohorts[name], condition, gA, gB)
        allA += list(A); allB += list(B)
        n, ngt, p, fold = paired_stat(A, B); statrows.append((name, n, ngt, p, fold))
        ax = axes[0, i]; c = PALETTE[i % len(PALETTE)]; fl = 1e-3
        for k in range(n):
            ax.plot([0, 1], [max(A[k], fl), max(B[k], fl)], color="#999", lw=0.8, zorder=1)
        ax.scatter(np.zeros(n), np.clip(A, fl, None), s=34, c=c, edgecolor="white", lw=0.5, zorder=3)
        ax.scatter(np.ones(n), np.clip(B, fl, None), s=34, c="#9aa0a6", edgecolor="white", lw=0.5, zorder=3)
        ax.set_yscale("log"); ax.set_xlim(-0.4, 1.4); ax.set_xticks([0, 1]); ax.set_xticklabels([gA, gB])
        ax.set_ylabel("%s TPM (log)" % gene); ax.set_title("%s (%s n=%d)" % (name, condition, n), fontsize=9.5)
        ax.text(0.5, 1.16, "%d/%d  %.0fx  P=%.3f" % (ngt, n, fold, p), transform=ax.transAxes,
                ha="center", va="top", fontsize=8, color="#333")
        ax.spines[["top", "right"]].set_visible(False)
    cn, cgt, cp, cfold = paired_stat(np.array(allA), np.array(allB))
    fig.suptitle("%s: %s vs %s - combined %s n=%d, %d/%d, P=%.4g"
                 % (gene, gA, gB, condition, cn, cgt, cn, cp), fontsize=10.5, fontweight="bold", y=1.04)
    for ext in ("png", "pdf", "svg"):
        fig.savefig("%s.%s" % (out, ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    with open(out + "_stats.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["cohort", "n", "%s>%s" % (gA, gB), "median_fold", "paired_wilcoxon_P"])
        for r in statrows:
            w.writerow([r[0], r[1], "%d/%d" % (r[2], r[1]), "%.2f" % r[4], "%.4g" % r[3]])
        w.writerow(["COMBINED", cn, "%d/%d" % (cgt, cn), "%.2f" % cfold, "%.4g" % cp])
    return {"per_cohort": statrows, "combined": (cn, cgt, cp, cfold)}
