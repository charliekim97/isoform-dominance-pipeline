#!/usr/bin/env python3
"""
04 - Contamination control (optional, config-driven).

Tests whether a target isoform group's signal could be explained by contamination from
another cell type. For each sample it computes marker scores from two gene panels and a
contamination score, then correlates the target-group abundance with that contamination
score (Spearman). A positive, significant correlation would suggest the target signal
tracks contamination; no positive correlation supports a genuine signal.

Requires: numpy, scipy, matplotlib.

config.json needs:
  "contamination_qc": {
    "target_group": "long",
    "marker_panels": {
      "tissue":      ["TTR","FOLR1","OTX2"],      # expected HIGH in the tissue of interest
      "contaminant": ["RBFOX3","SNAP25","GAD1"]    # expected HIGH only if contaminated
    }
  }

Inputs (repeat --markers / --target per cohort, names must match):
  --markers GSE228458=marker_tpm_GSE228458.csv   (CSV: donor + one TPM column per marker gene)
  --target  GSE228458=perdonor_GSE228458.csv     (from step 02; provides <target_group>_TPM)

Usage:
  python3 scripts/04_contamination_qc.py --config config.json \
    --markers GSE228458=markers_228.csv --target GSE228458=perdonor_228.csv \
    --markers GSE137619=markers_137.csv --target GSE137619=perdonor_137.csv \
    --out results/contamination
"""
import argparse, json, csv, os, math
import numpy as np
from scipy.stats import spearmanr
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

def log2p1(x): return math.log2(x + 1.0)

def load_markers(path, tissue, contaminant):
    out = {}
    for r in csv.DictReader(open(path)):
        t = np.mean([log2p1(float(r[g])) for g in tissue if g in r])
        c = np.mean([log2p1(float(r[g])) for g in contaminant if g in r])
        out[r["donor"]] = (t, c, (c / t if t > 0 else float("nan")))
    return out

def load_target(path, target_group):
    col = "%s_TPM" % target_group
    return {r["donor"]: float(r[col]) for r in csv.DictReader(open(path))}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--markers", action="append", required=True, help="NAME=marker_tpm.csv")
    ap.add_argument("--target", action="append", required=True, help="NAME=perdonor.csv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cfg = json.load(open(a.config))
    qc = cfg["contamination_qc"]
    tg = qc["target_group"]
    tissue = qc["marker_panels"]["tissue"]; contam = qc["marker_panels"]["contaminant"]
    markers = dict(kv.split("=", 1) for kv in a.markers)
    targets = dict(kv.split("=", 1) for kv in a.target)
    cohorts = list(markers)

    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 9, "pdf.fonttype": 42, "svg.fonttype": "none"})
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig, axes = plt.subplots(1, len(cohorts), figsize=(3.6 * len(cohorts), 3.4), squeeze=False)

    srows = []
    for i, name in enumerate(cohorts):
        m = load_markers(markers[name], tissue, contam)
        t = load_target(targets[name], tg)
        donors = [d for d in m if d in t]
        contam_score = np.array([m[d][1] for d in donors])
        ratio = np.array([m[d][2] for d in donors])
        target_val = np.array([t[d] for d in donors])
        rho, p = spearmanr(contam_score, target_val)
        srows.append((name, len(donors), rho, p, float(np.nanmedian(ratio))))
        ax = axes[0, i]; c = PALETTE[i % len(PALETTE)]
        ax.scatter(contam_score, target_val, s=36, c=c, edgecolor="white", lw=0.5, zorder=3)
        ax.set_xlabel("contamination score\n[mean log2(TPM+1), contaminant panel]")
        ax.set_ylabel("%s (%s) TPM" % (cfg.get("gene", "target"), tg))
        ax.set_title("%s (n=%d)" % (name, len(donors)), fontsize=9.5)
        ax.text(0.03, 0.97, "Spearman rho=%.2f\nP=%.2f" % (rho, p), transform=ax.transAxes,
                va="top", fontsize=8, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ccc", lw=0.6))
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("%s (%s) vs contamination score — no positive dependence = genuine signal"
                 % (cfg.get("gene", "target"), tg), fontsize=10, fontweight="bold", y=1.04)
    for ext in ("png", "pdf", "svg"):
        fig.savefig("%s.%s" % (a.out, ext), dpi=300, bbox_inches="tight")

    with open(a.out + "_scores.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["cohort", "n", "spearman_rho", "spearman_P", "median_contam_tissue_ratio"])
        for r in srows: w.writerow([r[0], r[1], "%.3f" % r[2], "%.4g" % r[3], "%.4f" % r[4]])

    print("=== contamination QC: %s (%s) vs contaminant score ===" % (cfg.get("gene", "target"), tg))
    for r in srows:
        print("  %-12s n=%d  rho=%+.3f  P=%.3f  contam/tissue ratio=%.3f" % (r[0], r[1], r[2], r[3], r[4]))
    print("wrote %s.{png,pdf,svg} + %s_scores.csv" % (a.out, a.out))

if __name__ == "__main__":
    main()
