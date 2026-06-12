#!/usr/bin/env python3
"""
03 — Paired stats + figure for isoform-group dominance (config-driven, generic).

Compares primary_comparison[0] vs [1] within each donor (paired Wilcoxon signed-rank),
per cohort and combined, on a chosen condition subset (e.g. control). Produces a
per-donor paired before-after figure (log y) and a stats table.

Requires: numpy, scipy, matplotlib.

Usage:
  python3 03_stats_and_figure.py \
      --config config.json \
      --condition control \
      --perdonor GSE228458=perdonor_GSE228458.csv \
      --perdonor GSE137619=perdonor_GSE137619.csv \
      --out results/LEPR_dominance
"""
import argparse, json, csv, os, numpy as np
from scipy.stats import wilcoxon
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt

PALETTE = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3', '#937860']

def load(path, condition, gA, gB):
    don, A, B = [], [], []
    for r in csv.DictReader(open(path)):
        if condition not in (None, '', 'all') and r['condition'] != condition:
            continue
        don.append(r['donor']); A.append(float(r['%s_TPM' % gA])); B.append(float(r['%s_TPM' % gB]))
    return don, np.array(A), np.array(B)

def stat(A, B):
    n = len(A); ngt = int(np.sum(A > B))
    try: p = wilcoxon(A, B, alternative='two-sided').pvalue
    except Exception: p = float('nan')
    with np.errstate(divide='ignore', invalid='ignore'):
        fold = np.median(A / B)
    return n, ngt, p, fold

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--condition', default='control')
    ap.add_argument('--perdonor', action='append', required=True, help='NAME=path.csv (repeatable)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    cfg = json.load(open(a.config))
    gene = cfg.get('gene', 'gene')
    gA, gB = cfg.get('primary_comparison', list(cfg['groups'])[:2])
    cohorts = {}
    for kv in a.perdonor:
        name, path = kv.split('=', 1); cohorts[name] = path

    mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 9, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)

    fig, axes = plt.subplots(1, len(cohorts), figsize=(3.6 * len(cohorts), 3.6), squeeze=False)
    allA, allB = [], []
    statrows = []
    for i, (name, path) in enumerate(cohorts.items()):
        don, A, B = load(path, a.condition, gA, gB)
        allA += list(A); allB += list(B)
        n, ngt, p, fold = stat(A, B); statrows.append((name, n, ngt, fold, p))
        ax = axes[0, i]; c = PALETTE[i % len(PALETTE)]; fl = 1e-3
        for k in range(n):
            ax.plot([0, 1], [max(A[k], fl), max(B[k], fl)], color='#999', lw=0.8, zorder=1)
        ax.scatter(np.zeros(n), np.clip(A, fl, None), s=34, c=c, edgecolor='white', lw=0.5, zorder=3)
        ax.scatter(np.ones(n), np.clip(B, fl, None), s=34, c='#9aa0a6', edgecolor='white', lw=0.5, zorder=3)
        ax.set_yscale('log'); ax.set_xlim(-0.4, 1.4); ax.set_xticks([0, 1]); ax.set_xticklabels([gA, gB])
        ax.set_ylabel('%s TPM (log)' % gene); ax.set_title('%s (%s n=%d)' % (name, a.condition, n), fontsize=9.5)
        ax.text(0.5, 1.16, '%d/%d  %.0fx  P=%.3f' % (ngt, n, fold, p), transform=ax.transAxes,
                ha='center', va='top', fontsize=8, color='#333')
        ax.spines[['top', 'right']].set_visible(False)
    # combined
    cn, cgt, cp, cfold = stat(np.array(allA), np.array(allB))
    fig.suptitle('%s: %s vs %s — combined %s n=%d, %d/%d, P=%.4g'
                 % (gene, gA, gB, a.condition, cn, cgt, cn, cp), fontsize=10.5, fontweight='bold', y=1.04)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig('%s.%s' % (a.out, ext), dpi=300, bbox_inches='tight')

    with open(a.out + '_stats.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['cohort', 'n', '%s>%s' % (gA, gB), 'median_fold', 'paired_wilcoxon_P'])
        for r in statrows: w.writerow([r[0], r[1], '%d/%d' % (r[2], r[1]), '%.2f' % r[3], '%.4g' % r[4]])
        w.writerow(['COMBINED', cn, '%d/%d' % (cgt, cn), '%.2f' % cfold, '%.4g' % cp])

    print('=== %s : %s vs %s (%s) ===' % (gene, gA, gB, a.condition))
    for r in statrows: print('  %-12s n=%d  %d/%d  fold=%.1fx  P=%.4g' % (r[0], r[1], r[2], r[1], r[3], r[4]))
    print('  COMBINED      n=%d  %d/%d  fold=%.1fx  P=%.4g' % (cn, cgt, cn, cfold, cp))
    print('wrote %s.{png,pdf,svg} + %s_stats.csv' % (a.out, a.out))

if __name__ == '__main__':
    main()
