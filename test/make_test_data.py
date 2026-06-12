#!/usr/bin/env python3
"""
Generate synthetic Salmon quant.sf files from real per-transcript LEPR TPM values,
so the whole pipeline can be validated end-to-end with NO external downloads.

The values below are the actual per-transcript TPMs from the published LEPR analysis
(control choroid plexus, GSE228458 + GSE137619). Running the pipeline on them must
reproduce: short>long in 5/5 and 6/6 donors; combined n=11, paired Wilcoxon P ~= 0.001.
"""
import os, csv

# donor: (ENST00000371060, ENST00000616738, ENST00000349533) TPM
DATA = {
    "GSE228458": [
        ("ctrl1", 6.521574, 36.556954, 1.314808),
        ("ctrl2", 1.250224,  8.637422, 0.349876),
        ("ctrl3", 2.238854, 13.496253, 0.704219),
        ("ctrl4", 3.684999, 15.015936, 0.721059),
        ("ctrl5", 0.495420,  6.599920, 1.440164),
    ],
    "GSE137619": [
        ("ctrl1", 9.130898,  9.898263, 0.351829),
        ("ctrl2", 5.896581,  5.928743, 0.216838),
        ("ctrl3", 7.516654, 20.317774, 1.033148),
        ("ctrl4", 1.804960,  3.059827, 0.143530),
        ("ctrl5", 3.730544,  2.848134, 0.203918),
        ("ctrl6", 2.582321,  4.337355, 0.166890),
    ],
}

def _write_quant(path, t1, t2, tl):
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, "quant.sf"), "w") as f:
        f.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        f.write("ENST00000371060.5\t3000\t2800\t%s\t100\n" % t1)
        f.write("ENST00000616738.1\t3000\t2800\t%s\t100\n" % t2)
        f.write("ENST00000349533.11\t4000\t3800\t%s\t10\n" % tl)

def generate(base):
    """Create quant/<cohort>/<donor>/quant.sf and sample_map_<cohort>.csv under `base`.
    Returns {cohort: {'quantdir':..., 'samplemap':...}}."""
    out = {}
    for cohort, rows in DATA.items():
        qd = os.path.join(base, "quant_%s" % cohort)
        for donor, t1, t2, tl in rows:
            _write_quant(os.path.join(qd, donor), t1, t2, tl)
        sm = os.path.join(base, "sample_map_%s.csv" % cohort)
        with open(sm, "w", newline="") as f:
            w = csv.writer(f); w.writerow(["donor", "condition", "SRR"])
            for donor, *_ in rows: w.writerow([donor, "control", "SYNTHETIC"])
        out[cohort] = {"quantdir": qd, "samplemap": sm}
    return out

TISSUE = ["TTR", "FOLR1", "OTX2", "AQP1"]          # expected HIGH (purity)
CONTAM = ["RBFOX3", "SNAP25", "MAP2", "GAD1"]      # expected LOW unless contaminated

def generate_markers(base):
    """Synthetic marker-gene TPM matrices: tissue markers high, contaminant markers low
    (mimics a pure tissue). Returns {cohort: path}."""
    out = {}
    for cohort, rows in DATA.items():
        path = os.path.join(base, "markers_%s.csv" % cohort)
        with open(path, "w", newline="") as f:
            w = csv.writer(f); w.writerow(["donor"] + TISSUE + CONTAM)
            for i, row in enumerate(rows):
                donor = row[0]
                tvals = [40000 + i * 3000, 1500 + i * 90, 800 + i * 40, 1200 + i * 70]
                cvals = [round(0.05 + 0.01 * i, 3), round(0.30 + 0.03 * i, 3),
                         round(0.10 + 0.02 * i, 3), round(0.02 + 0.01 * i, 3)]
                w.writerow([donor] + tvals + cvals)
        out[cohort] = path
    return out

if __name__ == "__main__":
    import sys
    base = sys.argv[1] if len(sys.argv) > 1 else "test/_generated"
    info = generate(base)
    mk = generate_markers(base)
    print("generated synthetic test data under", base)
    for c, p in info.items(): print(" ", c, p)
    for c, p in mk.items(): print(" markers", c, p)
