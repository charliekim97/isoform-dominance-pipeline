"""Download-free end-to-end self-test.

Generates synthetic Salmon quant.sf files from the per-transcript LEPR TPM values of the
choroid-plexus reanalysis, runs extract + stats (+ contamination), and asserts the reference
result is reproduced: combined n=11, 11/11 short>long, paired Wilcoxon P ~= 9.77e-4.
"""
import tempfile, os, csv, math, shutil
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


#: Order in which the cohort combinations are reported.
COMBINATIONS = ("pooled", "stouffer", "stratified_signed_rank")


def evaluate(base):
    """Run the pipeline on synthetic data under `base` and return the structured result.

    ``{"ok": bool, "checks": [...], "combinations": {...}, "headline_combination": str}``.
    Each check is ``{"name", "ok", "message", ...}`` plus the values it compared.  The
    three cohort combinations are reported, not checked: only the donor-pooled P is
    the reference result, and it is already a check.
    """
    info = generate(base)
    perdonor = {}
    for cohort, p in info.items():
        out = os.path.join(base, "perdonor_%s.csv" % cohort)
        extract.run(CONFIG, p["quantdir"], p["samplemap"], cohort, out)
        perdonor[cohort] = out
    res = stats.run(CONFIG, "control", perdonor, os.path.join(base, "result"))
    checks = []

    def check(name, good, message, **values):
        checks.append(dict(name=name, ok=bool(good),
                           message="[%s] %s" % ("OK" if good else "FAIL", message), **values))

    by = {d["cohort"]: d for d in res["detail"]}
    for coh in ("GSE228458", "GSE137619"):
        n, ngt = EXPECT[coh]
        d = by[coh]
        check(coh, d["n"] == n and d["n_greater"] == ngt,
              "%s %d/%d short>long  P=%.4g (%s)"
              % (coh, d["n_greater"], d["n"], d["p"], d["wilcoxon_method"]),
              n=d["n"], n_greater=d["n_greater"], expect=[n, ngt], p=d["p"],
              wilcoxon_method=d["wilcoxon_method"])
    cn, cgt, cp, _ = res["combined"]
    en, eng, ep = EXPECT["COMBINED"]
    check("COMBINED", cn == en and cgt == eng and abs(cp - ep) < 1e-4,
          "COMBINED %d/%d P=%.4g (expect %.4g)" % (cgt, cn, cp, ep),
          n=cn, n_greater=cgt, p=cp, expect=[en, eng, ep],
          wilcoxon_method=res["pooled"]["wilcoxon_method"])
    fig_ok = os.path.exists(os.path.join(base, "result.png"))
    check("figure", fig_ok, "dominance figure produced")
    rows = contamination.run(CONFIG, {c: info[c]["markers"] for c in info},
                             {c: perdonor[c] for c in info}, os.path.join(base, "qc"))
    ratios = [r[4] for r in rows]
    qc_ok = os.path.exists(os.path.join(base, "qc.png")) and all(x < 0.2 for x in ratios)
    check("contamination_qc", qc_ok, "contamination-QC purity ratios %s < 0.2"
          % [round(x, 3) for x in ratios], ratios=ratios, threshold=0.2)

    combinations = {}
    for key in COMBINATIONS:
        c = res["combination"][key]
        combinations[key] = {"label": stats.COMBINATION_LABELS[key], "k": c["k"],
                             "p": c["p"],
                             "z": c["z"] if math.isfinite(c["z"]) else None}
    combinations["pooled"]["wilcoxon_method"] = res["pooled"]["wilcoxon_method"]
    return {"ok": all(c["ok"] for c in checks), "checks": checks,
            "combinations": combinations,
            "headline_combination": res["headline_combination"]}


def run(base):
    """Run the pipeline on synthetic data under `base`. Returns (ok, messages).

    The v2.1 shape, kept for callers; :func:`evaluate` returns the structure.
    """
    r = evaluate(base)
    return r["ok"], [c["message"] for c in r["checks"]]


def result():
    """Run the self-test in a temporary directory and return :func:`evaluate`'s dict."""
    work = tempfile.mkdtemp(prefix="idp_selftest_")
    try:
        return evaluate(work)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    r = result()
    for c in r["checks"]:
        print("  " + c["message"])
    print("\n  cohort combinations (reported; only donor-pooled is checked above):")
    for key in COMBINATIONS:
        c = r["combinations"][key]
        line = "    %-24s k=%d  P=%.4g" % (c["label"], c["k"], c["p"])
        if key == "pooled":
            line += "  (paired Wilcoxon, %s)" % c["wilcoxon_method"]
        if key == r["headline_combination"]:
            line += "  <- headline"
        print(line)
    print("\n%s" % ("PASS - reproduces the reference LEPR result." if r["ok"] else "FAIL"))
    return 0 if r["ok"] else 1
