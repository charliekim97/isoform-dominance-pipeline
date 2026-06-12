#!/usr/bin/env python3
"""
End-to-end self-test (no downloads). Generates synthetic quant.sf from real LEPR
per-transcript TPM, runs the extraction and stats steps, and asserts the published
result is reproduced: combined n=11, 11/11 short>long, paired Wilcoxon P ~= 0.000977.

Run from the repository root:
    python3 test/run_test.py
Exit code 0 = PASS.
"""
import os, sys, csv, tempfile, subprocess, shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "test"))
import make_test_data  # noqa: E402

EXPECT = {  # cohort: (n, n_gt, P)  and combined
    "GSE228458": (5, 5, 0.0625),
    "GSE137619": (6, 6, 0.03125),
    "COMBINED":  (11, 11, 0.0009766),
}
PTOL = 1e-3

def main():
    work = tempfile.mkdtemp(prefix="idp_test_")
    try:
        info = make_test_data.generate(work)
        cfg = os.path.join(ROOT, "config.example.json")
        perdonor_args = []
        for cohort, p in info.items():
            out = os.path.join(work, "perdonor_%s.csv" % cohort)
            subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "02_extract_isoform_tpm.py"),
                            "--config", cfg, "--quantdir", p["quantdir"],
                            "--samplemap", p["samplemap"], "--cohort", cohort, "--out", out],
                           check=True)
            perdonor_args += ["--perdonor", "%s=%s" % (cohort, out)]
        prefix = os.path.join(work, "result")
        subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "03_stats_and_figure.py"),
                        "--config", cfg, "--condition", "control", "--out", prefix] + perdonor_args,
                       check=True)

        got = {}
        for r in csv.DictReader(open(prefix + "_stats.csv")):
            n = int(r["n"]); ngt = int(r[[k for k in r if ">" in k][0]].split("/")[0])
            got[r["cohort"]] = (n, ngt, float(r["paired_wilcoxon_P"]))

        ok = True
        for k, (n, ngt, p) in EXPECT.items():
            gn, gng, gp = got.get(k, (None, None, None))
            good = (gn == n and gng == ngt and abs(gp - p) < PTOL)
            print("  [%s] %-9s expected n=%d %d/%d P=%.4g | got n=%s %s/%s P=%.4g"
                  % ("OK" if good else "FAIL", k, n, ngt, n, p, gn, gng, gn, gp))
            ok = ok and good

        fig_ok = os.path.exists(prefix + ".png")
        print("  [%s] figure produced" % ("OK" if fig_ok else "FAIL"))
        ok = ok and fig_ok

        # --- step 04: contamination-QC module ---
        mk = make_test_data.generate_markers(work)
        qc_args = []
        for cohort in info:
            qc_args += ["--markers", "%s=%s" % (cohort, mk[cohort]),
                        "--target", "%s=%s" % (cohort, os.path.join(work, "perdonor_%s.csv" % cohort))]
        qc_prefix = os.path.join(work, "qc")
        subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "04_contamination_qc.py"),
                        "--config", cfg, "--out", qc_prefix] + qc_args, check=True)
        ratios = [float(r["median_contam_tissue_ratio"]) for r in csv.DictReader(open(qc_prefix + "_scores.csv"))]
        qc_ok = os.path.exists(qc_prefix + ".png") and all(x < 0.2 for x in ratios)
        print("  [%s] contamination-QC module runs; purity ratios %s (< 0.2)"
              % ("OK" if qc_ok else "FAIL", [round(x, 3) for x in ratios]))
        ok = ok and qc_ok

        print("\n%s" % ("PASS — pipeline reproduces the published LEPR result." if ok else "FAIL"))
        return 0 if ok else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)

if __name__ == "__main__":
    sys.exit(main())
