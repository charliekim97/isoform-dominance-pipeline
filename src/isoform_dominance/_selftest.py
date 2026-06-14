"""Download-free end-to-end self-test.

Generates synthetic Salmon quant.sf files from the real per-transcript LEPR TPM values of
the published analysis, runs extract + stats (+ contamination), and asserts the published
result is reproduced: combined n=11, 11/11 short>long, paired Wilcoxon P ~= 9.77e-4.
"""
import tempfile, os, csv, shutil
from . import extract, stats, contamination

# donor: (ENST00000371060, ENST00000616738, ENST00000349533) TPM  (control choroid plexus)
DATA = {
    "GSE228458": [("ctrl1", 6.521574, 36.556954, 1.314808), ("ctrl2", 1.250224, 8.637422, 0.349876),
                  ("ctrl3", 2.238854, 13.496253, 0.704219), ("ctrl4", 3.684999, 15.015936, 0.721059),
                  ("ctrl5", 0.495420, 6.599920, 1.440164)],
    "GSE137619": [("ctrl1", 9.130898, 9.898263, 0.351829), ("ctrl2", 5.896581, 5.928743, 0.216838),
                  ("ctrl3", 7.516654, 20.317774, 1.033148), ("ctrl4", 1.804960, 3.059827, 0.143530),
                  ("ctrl5", 3.730544, 2.848134, 0.203918), ("ctrl6", 2.582321, 4.337355, 0.166890)],
}
CONFIG = {
    "gene": "LEPR",
    "groups": {"short": ["ENST00000371060", "ENST00000616738"], "long": ["ENST00000349533"]},
    "primary_comparison": ["short", "long"],
    "contamination_qc": {"target_group": "long",
                         "marker_panels": {"tissue": ["TTR", "FOLR1", "OTX2", "AQP1"],
                                           "contaminant": ["RBFOX3", "SNAP25", "MAP2", "GAD1"]}},
}
EXPECT = {"GSE228458": (5, 5), "GSE137619": (6, 6), "COMBINED": (11, 11, 9.7656e-4)}


def _write_quant(path, t1, t2, tl):
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, "quant.sf"), "w") as f:
        f.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        f.write("ENST00000371060.5\t3000\t2800\t%s\t100\n" % t1)
        f.write("ENST00000616738.1\t3000\t2800\t%s\t100\n" % t2)
        f.write("ENST00000349533.11\t4000\t3800\t%s\t10\n" % tl)


def generate(base):
    out = {}
    for cohort, rows in DATA.items():
        qd = os.path.join(base, "quant_%s" % cohort)
        for donor, t1, t2, tl in rows:
            _write_quant(os.path.join(qd, donor), t1, t2, tl)
        sm = os.path.join(base, "sample_map_%s.csv" % cohort)
        with open(sm, "w", newline="") as f:
            w = csv.writer(f); w.writerow(["donor", "condition", "SRR"])
            for donor, *_ in rows:
                w.writerow([donor, "control", "SYNTHETIC"])
        mk = os.path.join(base, "markers_%s.csv" % cohort)
        with open(mk, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["donor", "TTR", "FOLR1", "OTX2", "AQP1", "RBFOX3", "SNAP25", "MAP2", "GAD1"])
            for i, (donor, *_rest) in enumerate(rows):
                w.writerow([donor, 40000 + i * 3000, 1500 + i * 90, 800 + i * 40, 1200 + i * 70,
                            round(0.05 + 0.01 * i, 3), round(0.30 + 0.03 * i, 3),
                            round(0.10 + 0.02 * i, 3), round(0.02 + 0.01 * i, 3)])
        out[cohort] = {"quantdir": qd, "samplemap": sm, "markers": mk}
    return out


def run(base):
    """Run the pipeline on synthetic data under `base`. Returns (ok, messages)."""
    info = generate(base)
    perdonor = {}
    for cohort, p in info.items():
        out = os.path.join(base, "perdonor_%s.csv" % cohort)
        extract.run(CONFIG, p["quantdir"], p["samplemap"], cohort, out)
        perdonor[cohort] = out
    res = stats.run(CONFIG, "control", perdonor, os.path.join(base, "result"))
    msgs, ok = [], True
    by = {r[0]: r for r in res["per_cohort"]}
    for coh, (n, ngt) in [("GSE228458", EXPECT["GSE228458"]), ("GSE137619", EXPECT["GSE137619"])]:
        gn, gng = by[coh][1], by[coh][2]
        good = gn == n and gng == ngt; ok &= good
        msgs.append("[%s] %s %d/%d short>long" % ("OK" if good else "FAIL", coh, gng, gn))
    cn, cgt, cp, _ = res["combined"]
    en, eng, ep = EXPECT["COMBINED"]
    cgood = cn == en and cgt == eng and abs(cp - ep) < 1e-4; ok &= cgood
    msgs.append("[%s] COMBINED %d/%d P=%.4g (expect %.4g)" % ("OK" if cgood else "FAIL", cgt, cn, cp, ep))
    fig_ok = os.path.exists(os.path.join(base, "result.png")); ok &= fig_ok
    msgs.append("[%s] dominance figure produced" % ("OK" if fig_ok else "FAIL"))
    rows = contamination.run(CONFIG, {c: info[c]["markers"] for c in info},
                             {c: perdonor[c] for c in info}, os.path.join(base, "qc"))
    ratios = [r[4] for r in rows]
    qc_ok = os.path.exists(os.path.join(base, "qc.png")) and all(x < 0.2 for x in ratios); ok &= qc_ok
    msgs.append("[%s] contamination-QC purity ratios %s < 0.2"
                % ("OK" if qc_ok else "FAIL", [round(x, 3) for x in ratios]))
    return ok, msgs


def main():
    work = tempfile.mkdtemp(prefix="idp_selftest_")
    try:
        ok, msgs = run(work)
        for m in msgs:
            print("  " + m)
        print("\n%s" % ("PASS - reproduces the published LEPR result." if ok else "FAIL"))
        return 0 if ok else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)
