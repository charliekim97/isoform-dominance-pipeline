"""Grade a TPM run against the 3.0 gate of ``validation/README.md``.

    python validation/gate_3_0.py c1_d1_tpm.json

The run is ``run_simulation.py run`` at the stored design and uniform coverage, quantified by
Salmon 2.7.0 with no option beyond those the harness passes; its ``genes_tpm`` rows are the
class ratio of summed TPM against the ratio summed TPM would read had every fragment gone to
the transcript it was drawn from.  The predictions are the ones fixed with the gate,
``gate_3_0_predictions.json``, whose SHA-256 the gate's text gives; nothing here recomputes
them.  This script applies the gate's steps 3 to 6 as written, in its order:

- the primary set first: every one of its 7 genes with a |bias| below its predicted SE, at
  least 6 of the 7 with a replicate SD below it, and the stop rule not triggered;
- the stop rule: SMN2 or VEGFA with a |bias| or an SD at or above its predicted SE;
- the control only when the primary set passes: of its 39 genes, a |bias| below the predicted
  SE for at least 38 and an SD below it for at least 31.

A gene without a row (fewer than 3 replicates that count) fails both criteria.  The exit
status is 0 when the gate passes (3.0 makes ``abundance`` the default), 1 when it does not,
and 2 when the run or the predictions are not the ones the gate names; then nothing is graded.
A run is the one the gate names when its design, quantifier and Salmon version are the stored
run's, its libraries held every gene graded, and its ``NumReads`` rows reproduce the stored
uniform run, ``data/sim/c1_d1.json``, within :data:`REPRODUCE_TOLERANCE`: the same seeds,
reads, Salmon and index give the same rows, so a larger difference is a run that was not
made as stored, which is found here without reading a TPM row.
"""
import argparse
import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PREDICTIONS = os.path.join(HERE, "gate_3_0_predictions.json")
README = os.path.join(HERE, "README.md")

#: the gate's numbers, as validation/README.md states them (tests/test_gate_3_0.py holds the
#: text to them)
STOP_GENES = ("SMN2", "VEGFA")
PRIMARY_SD_AT_LEAST = 6
CONTROL_BIAS_AT_LEAST = 38
CONTROL_SD_AT_LEAST = 31
#: the design the gate repeats: the stored one, at uniform coverage
DESIGN = {"depth_scale": 1.0, "ratio": None, "tpm": 10.0, "bias_5p": 0.0,
          "mean_efflen": 1500.0, "replicates": 30, "salmon_flags": "",
          "quantifier": "salmon"}
SALMON_VERSION = "2.7.0"
#: the stored uniform-coverage run the NumReads rows must reproduce
STORED_RUN = os.path.join(HERE, "data", "sim", "c1_d1.json")
#: the largest difference, in log2, of a NumReads row's bias, SD or RMSE from the stored
#: run's: a rerun on the machine that made it gave identical rows (``chk_d1.json``), and this
#: leaves room only for a quantifier's last digits
REPRODUCE_TOLERANCE = 0.01


class NotGradeable(Exception):
    """The run or the predictions are not the ones the gate names."""


def frozen_predictions(path=PREDICTIONS, readme=README):
    """The fixed predictions, after checking their bytes against the SHA-256 the gate's
    text gives."""
    with open(path, "rb") as f:
        data = f.read()
    with open(readme, encoding="utf-8") as f:
        text = f.read()
    m = re.search(r"`gate_3_0_predictions\.json` beside this file.*?SHA-256 `([0-9a-f]{64})`",
                  text, re.S)
    if not m:
        raise NotGradeable("%s gives no SHA-256 for gate_3_0_predictions.json" % readme)
    got = hashlib.sha256(data).hexdigest()
    if got != m.group(1):
        raise NotGradeable("gate_3_0_predictions.json has SHA-256 %s, not the %s the gate "
                           "names" % (got, m.group(1)))
    return json.loads(data.decode("utf-8"))


def check_run(run, predictions):
    """Refuse a run that is not the one the gate names."""
    design = run.get("design") or {}
    for key, want in DESIGN.items():
        if design.get(key) != want:
            raise NotGradeable("the run's design has %s %r, not %r" % (key, design.get(key),
                                                                       want))
    version = (run.get("salmon") or {}).get("version")
    if version != SALMON_VERSION:
        raise NotGradeable("the run was quantified by Salmon %r, not %s"
                           % (version, SALMON_VERSION))
    if "genes_tpm" not in run:
        raise NotGradeable("the run has no genes_tpm rows: it was made by a harness that "
                           "does not record summed TPM")
    wanted = set(predictions["primary"]) | set(predictions["control"])
    simulated = set(design.get("genes_simulated") or ())
    if not wanted <= simulated:
        raise NotGradeable("the run's libraries did not hold %s"
                           % ", ".join(sorted(wanted - simulated)))


def check_reproduces(run, stored=STORED_RUN):
    """The largest difference of the run's NumReads rows from the stored run's; refuses a
    run whose rows are for other genes, count other replicates, or differ by more than
    :data:`REPRODUCE_TOLERANCE`."""
    with open(stored, encoding="utf-8") as f:
        want = {r["gene"]: r for r in json.load(f)["genes"]}
    got = {r["gene"]: r for r in run.get("genes") or ()}
    if set(got) != set(want):
        raise NotGradeable("the run's NumReads rows are for %d genes, the stored run's for %d "
                           "(%s)" % (len(got), len(want),
                                     ", ".join(sorted(set(got) ^ set(want)))))
    worst = 0.0
    for gene, w in sorted(want.items()):
        r = got[gene]
        if r["n"] != w["n"]:
            raise NotGradeable("%s's NumReads row counts %d replicates, the stored run's %d"
                               % (gene, r["n"], w["n"]))
        for key in ("bias", "realised_sd", "rmse"):
            worst = max(worst, abs(r[key] - w[key]))
    if worst > REPRODUCE_TOLERANCE:
        raise NotGradeable("the run's NumReads rows differ from the stored run's by up to "
                           "%.4f log2, more than %.2f: it was not made as stored"
                           % (worst, REPRODUCE_TOLERANCE))
    return worst


def _grade(rows, predicted):
    """Per gene ``(row or None, bias ok, SD ok)``; a gene without a row fails both."""
    out = {}
    for gene, se in predicted.items():
        r = rows.get(gene)
        if r is None:
            out[gene] = (None, False, False)
        else:
            out[gene] = (r, abs(r["bias"]) < se, r["realised_sd"] < se)
    return out


def grade(run, predictions):
    """The gate's verdict on ``run``: a dict with the per-gene grades, ``stopped``,
    ``primary_passes``, the control's counts when it is graded, and ``passes``."""
    rows = {r["gene"]: r for r in run["genes_tpm"]}
    primary = _grade(rows, predictions["primary"])
    stopped = [g for g in STOP_GENES if not (primary[g][1] and primary[g][2])]
    sd_ok = sum(1 for _, _, s in primary.values() if s)
    primary_passes = (all(b for _, b, _ in primary.values())
                      and sd_ok >= PRIMARY_SD_AT_LEAST and not stopped)
    out = {"primary": primary, "primary_sd_ok": sd_ok, "stopped": stopped,
           "primary_passes": primary_passes, "control": None, "passes": False}
    if primary_passes:
        control = _grade(rows, predictions["control"])
        n_bias = sum(1 for _, b, _ in control.values() if b)
        n_sd = sum(1 for _, _, s in control.values() if s)
        out["control"] = {"grades": control, "bias_ok": n_bias, "sd_ok": n_sd,
                          "passes": (n_bias >= CONTROL_BIAS_AT_LEAST
                                     and n_sd >= CONTROL_SD_AT_LEAST)}
        out["passes"] = out["control"]["passes"]
    return out


def _line(gene, se, grade_):
    r, b, s = grade_
    if r is None:
        return "  %-9s predicted SE %.4f   no row (fewer than 3 replicates that count)" % (gene,
                                                                                         se)
    return ("  %-9s predicted SE %.4f   |bias| %.4f %s   SD %.4f %s   n %d"
            % (gene, se, abs(r["bias"]), "below" if b else "AT OR ABOVE",
               r["realised_sd"], "below" if s else "AT OR ABOVE", r["n"]))


def report(result, predictions):
    print("The 3.0 gate (validation/README.md), steps 3 to 6.")
    print("\nThe primary set (step 3):")
    for gene, se in predictions["primary"].items():
        print(_line(gene, se, result["primary"][gene]))
    print("  SD below the predicted SE for %d of %d (at least %d needed)"
          % (result["primary_sd_ok"], len(predictions["primary"]), PRIMARY_SD_AT_LEAST))
    if result["stopped"]:
        print("\nThe stop rule (step 4) stops the change: %s. Why a quantifier's TPM departs "
              "from an estimable abundance is found first." % ", ".join(result["stopped"]))
    if not result["primary_passes"]:
        print("\nThe primary set fails; the control is not graded.")
    else:
        c = result["control"]
        print("\nThe primary set passes. The control (step 6):")
        for gene, se in predictions["control"].items():
            print(_line(gene, se, c["grades"][gene]))
        print("  |bias| below the predicted SE for %d of %d (at least %d needed); SD below it "
              "for %d (at least %d)" % (c["bias_ok"], len(predictions["control"]),
                                        CONTROL_BIAS_AT_LEAST, c["sd_ok"],
                                        CONTROL_SD_AT_LEAST))
    print("\n%s" % ("GATE PASSES: 3.0 makes `abundance` the default." if result["passes"]
                    else "GATE FAILS: 3.0 keeps `fragment` as the default."))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="the run_simulation.py run output to grade")
    ap.add_argument("--predictions", default=PREDICTIONS)
    ap.add_argument("--readme", default=README)
    a = ap.parse_args(argv)
    try:
        predictions = frozen_predictions(a.predictions, a.readme)
        with open(a.run, encoding="utf-8") as f:
            run = json.load(f)
        check_run(run, predictions)
        worst = check_reproduces(run)
    except NotGradeable as e:
        print("NOT GRADED: %s" % e, file=sys.stderr)
        return 2
    print("The run's NumReads rows reproduce the stored uniform run (largest difference %.2g "
          "log2).\n" % worst)
    result = grade(run, predictions)
    report(result, predictions)
    return 0 if result["passes"] else 1


if __name__ == "__main__":
    sys.exit(main())
