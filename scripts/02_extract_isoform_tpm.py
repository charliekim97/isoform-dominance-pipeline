#!/usr/bin/env python3
"""
02 — Extract per-donor isoform-GROUP TPM from Salmon quant.sf (config-driven, generic).

Define any transcript groups in config.json (e.g. short vs long, variant A vs B).
Standard-library only -> runs on an HPC login node without numpy/scipy.

Usage:
  python3 02_extract_isoform_tpm.py \
      --config config.json \
      --quantdir /path/to/quant   (contains <donor>/quant.sf) \
      --samplemap sample_map.csv  (columns: donor,condition[,SRR]) \
      --cohort GSE228458 \
      --out perdonor_GSE228458.csv

config.json (example):
{
  "gene": "LEPR",
  "groups": { "short": ["ENST00000371060","ENST00000616738"],
              "long":  ["ENST00000349533"] },
  "primary_comparison": ["short","long"]
}
"""
import argparse, json, csv, os, glob

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--quantdir', required=True)
    ap.add_argument('--samplemap', required=True)
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    cfg = json.load(open(a.config))
    groups = cfg['groups']                       # {group_name: [ENST, ...]}
    prim = cfg.get('primary_comparison', list(groups)[:2])
    tx2grp = {}                                  # transcript(no version) -> group
    for g, txs in groups.items():
        for t in txs:
            tx2grp[t.split('.')[0]] = g

    cond = {}
    for r in csv.DictReader(open(a.samplemap)):
        cond[r['donor']] = r.get('condition', 'NA')

    rows = []
    for q in sorted(glob.glob(os.path.join(a.quantdir, '*', 'quant.sf'))):
        donor = os.path.basename(os.path.dirname(q))
        gt = {g: 0.0 for g in groups}
        for r in csv.DictReader(open(q), delimiter='\t'):
            tid = r['Name'].split('.')[0]
            g = tx2grp.get(tid)
            if g:
                gt[g] += float(r['TPM'])
        rows.append((donor, cond.get(donor, 'NA'), gt))

    with open(a.out, 'w', newline='') as f:
        w = csv.writer(f)
        head = ['cohort', 'donor', 'condition'] + ['%s_TPM' % g for g in groups]
        if len(prim) == 2:
            head.append('%s_fraction' % prim[0])
        w.writerow(head)
        for donor, c, gt in rows:
            row = [a.cohort, donor, c] + ['%.4f' % gt[g] for g in groups]
            if len(prim) == 2:
                s = gt[prim[0]] + gt[prim[1]]
                row.append('%.4f' % (gt[prim[0]] / s) if s > 0 else 'NA')
            w.writerow(row)
    print('wrote', a.out, 'n=', len(rows))

if __name__ == '__main__':
    main()
