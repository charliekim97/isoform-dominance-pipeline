"""Freeze, run and grade the three prospective G2 conditions.

Use `python -m validation.gate_g2` from the checkout. Freeze reads the reference
and computes predictions before a library exists. Run and grade take the SHA-256
recorded at freeze time. Grade exits 0 for pass, 1 for failure, and 2 for a
design or provenance mismatch.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

from isoform_dominance import io

_HARNESS_SPEC = importlib.util.spec_from_file_location(
    "_g2_harness", Path(__file__).with_name("run_simulation.py"))
harness = importlib.util.module_from_spec(_HARNESS_SPEC)
_HARNESS_SPEC.loader.exec_module(harness)

FORMAT = "isoform-dominance-g2-v1"
INPUTS_SHA256 = "ebbb3654cb1ec816eb71e561aa2481f681a15654db339eb6f4614fc7d565455b"
REFERENCE_TEXT_SHA256 = "50993c71f1ef9a5b27e2e492bcfba032721c3e784c545a74be6ba02aefa1720a"
CONDITIONS = {"flat": None, "a4": 4.0, "b4": 0.25}
DESIGN = {"depth_scale": 1.0, "ratio": None, "tpm": 10.0, "bias_5p": 0.0,
          "mean_efflen": 1500.0, "replicates": 30, "salmon_flags": "",
          "quantifier": "salmon"}
SALMON_VERSION = "2.7.0"
PRIMARY = ["CANX", "GNAS", "NTRK2", "NTRK3", "TNFRSF1B", "VEGFA"]
REPORT_ONLY = ["SMN2", "TSC1"]
RULES = {"primary_bias": 6, "primary_sd": 5, "control_bias": 37, "control_sd": 32,
         "gross_bias_se": 2.0, "minimum_usable_replicates": 3}


def _text_sha256(path):
    digest = hashlib.sha256()
    with io.open_text(path) as source:
        for chunk in iter(lambda: source.read(2 ** 20), ""):
            digest.update(chunk.encode("utf-8"))
    return digest.hexdigest()


def _check_inputs(inputs, reference):
    if io.file_sha256(inputs) != INPUTS_SHA256:
        raise ValueError("G2 requires the archived 49-gene sweep bundle")
    if _text_sha256(reference) != REFERENCE_TEXT_SHA256:
        raise ValueError("G2 requires the exact Ensembl 116 reference-chromosome FASTA")


def _seed_offset(condition):
    return 1_000_003 * (list(CONDITIONS).index(condition) + 1)


def _seeds(condition):
    return [7919 * rep + 13 + _seed_offset(condition) for rep in range(DESIGN["replicates"])]


def _has_se(row):
    se = row["predicted_se"]
    return se is not None and math.isfinite(se) and se > 0


def _groups(genes, predicted):
    abundance = {g for g, row in predicted["abundance"].items() if _has_se(row)}
    fragment = {g for g, row in predicted["fragment"].items() if _has_se(row)}
    return {"primary": sorted(abundance - fragment), "control": sorted(abundance & fragment),
            "report_only": sorted(set(genes) - abundance)}


def _execution_identity():
    sources = {
        "gate": __file__, "harness": harness.__file__, "simulator": harness.SIM,
        "model": harness.I.__file__, "index_scope": harness.I.index_scope.__file__,
        "io": io.__file__,
    }
    return {"numpy_version": harness.np.__version__,
            "code_sha256": {name: io.file_sha256(path) for name, path in sources.items()}}


def _plan(genes, competitor_ids, predictions):
    plan = {"format": FORMAT, "inputs_sha256": INPUTS_SHA256,
            "reference_text_sha256": REFERENCE_TEXT_SHA256, "salmon_version": SALMON_VERSION,
            "design": DESIGN, "rules": RULES, "execution": _execution_identity(),
            "genes": sorted(genes),
            "competitor_transcripts": sorted(competitor_ids),
            "groups": _groups(genes, predictions["flat"]),
            "conditions": {name: {"abundance_ratio": ratio, "seed_offset": _seed_offset(name),
                                  "seeds": _seeds(name), "predictions": predictions[name]}
                           for name, ratio in CONDITIONS.items()}}
    _check_plan(plan)
    return plan


def _check_plan(plan):
    for key, expected in (("format", FORMAT), ("inputs_sha256", INPUTS_SHA256),
                          ("reference_text_sha256", REFERENCE_TEXT_SHA256),
                          ("salmon_version", SALMON_VERSION), ("design", DESIGN),
                          ("rules", RULES)):
        if plan.get(key) != expected:
            raise ValueError(f"the G2 plan has different {key}")
    genes, competitors = plan["genes"], plan["competitor_transcripts"]
    if len(genes) != 49 or genes != sorted(set(genes)):
        raise ValueError("the G2 plan must contain 49 distinct genes")
    if len(competitors) != 3793 or competitors != sorted(set(competitors)):
        raise ValueError("the G2 plan must contain the 3,793 additional transcript origins")
    groups = plan["groups"]
    if groups["primary"] != PRIMARY or groups["report_only"] != REPORT_ONLY or len(groups["control"]) != 41:
        raise ValueError(f"G2 eligibility differs from the proposed 6/41/2 panel: {groups}")
    if set(plan["conditions"]) != set(CONDITIONS):
        raise ValueError("the G2 plan must contain all three conditions")
    for name, ratio in CONDITIONS.items():
        condition = plan["conditions"][name]
        if (condition["abundance_ratio"] != ratio or condition["seed_offset"] != _seed_offset(name)
                or condition["seeds"] != _seeds(name)):
            raise ValueError(f"the {name} condition has a different mix or seeds")
        predictions = condition["predictions"]
        if (set(predictions) != {"abundance", "fragment"}
                or any(set(rows) != set(genes) for rows in predictions.values())
                or _groups(genes, predictions) != groups):
            raise ValueError(f"the {name} condition changes the panel or its eligibility")


def freeze(inputs, reference):
    _check_inputs(inputs, reference)
    source = harness._sim().read_inputs(inputs=str(inputs))
    backgrounds, competitors = harness._reference_background(source, str(reference))
    predictions = {
        name: harness._predictions_at_abundance_ratio(
            source, DESIGN["tpm"], DESIGN["depth_scale"], DESIGN["mean_efflen"], ratio,
            backgrounds=backgrounds, competitors=competitors)
        for name, ratio in CONDITIONS.items()
    }
    return _plan(source, competitors, io.finite_json(predictions))


def _read_frozen(path, expected_sha):
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != expected_sha:
        raise ValueError("prediction file SHA-256 differs from the frozen digest")
    plan = json.loads(content)
    _check_plan(plan)
    return plan


def _provenance(sha, name):
    return {"plan_sha256": sha, "condition": name, "reference_text_sha256": REFERENCE_TEXT_SHA256}


def run_condition(plan, sha, name, inputs, reference, index, work, out, threads):
    if work.exists() or out.exists():
        raise ValueError("use new work and result paths for each G2 condition")
    if not out.parent.is_dir() or out.resolve().is_relative_to(work.resolve()):
        raise ValueError("the result needs an existing parent outside the replicate work directory")
    if plan["execution"] != _execution_identity():
        raise ValueError("the running code or NumPy version differs from the frozen plan")
    identity = harness._reference_identity(index, reference)
    _check_inputs(inputs, reference)
    if identity["index"]["salmon_version"] != SALMON_VERSION:
        raise ValueError(f"G2 needs an index from Salmon {SALMON_VERSION}")
    version = subprocess.run(["salmon", "--version"], check=True, text=True,
                             capture_output=True).stdout.strip().split()[-1]
    if version != SALMON_VERSION:
        raise ValueError(f"G2 needs Salmon {SALMON_VERSION}, found {version}")
    case = plan["conditions"][name]
    argv = ["run", "--inputs", str(inputs), "--background-fasta", str(reference),
            "--index", str(index), "--work", str(work), "--out", str(out),
            "--threads", str(threads), "--seed-offset", str(case["seed_offset"]),
            "--depth-scale", str(DESIGN["depth_scale"]), "--tpm", str(DESIGN["tpm"]),
            "--mean-efflen", str(DESIGN["mean_efflen"]), "--replicates", str(DESIGN["replicates"])]
    if case["abundance_ratio"] is not None:
        argv += ["--abundance-ratio", str(case["abundance_ratio"])]
    work.mkdir(parents=True, exist_ok=False)
    harness.main(argv, frozen={"predictions": case["predictions"],
                               "verified_fasta_sha256": identity["fasta_sha256"],
                               "competitor_transcripts": plan["competitor_transcripts"],
                               "provenance": _provenance(sha, name)})


def _grade_rows(plan, name, rows):
    case = plan["conditions"][name]
    se = case["predictions"]["abundance"]
    results = {}
    for group in ("primary", "control"):
        by_gene, gross = {}, []
        for gene in plan["groups"][group]:
            row = rows.get(gene)
            limit = se[gene]["predicted_se"]
            evaluable = row is not None and row["n"] >= RULES["minimum_usable_replicates"]
            by_gene[gene] = {"bias_below_se": evaluable and abs(row["bias"]) < limit,
                             "sd_below_se": evaluable and row["realised_sd"] < limit,
                             "beyond_linear": se[gene]["beyond_linear"],
                             "usable_replicates": row["n"] if row else 0}
            if evaluable and abs(row["bias"]) >= RULES["gross_bias_se"] * limit:
                gross.append(gene)
        bias = sum(r["bias_below_se"] for r in by_gene.values())
        sd = sum(r["sd_below_se"] for r in by_gene.values())
        results[group] = {"bias_below_se": bias, "sd_below_se": sd, "gross": gross,
                          "passes": bias >= RULES[f"{group}_bias"]
                                    and sd >= RULES[f"{group}_sd"] and not gross,
                          "genes": by_gene}
    return {"passes": all(r["passes"] for r in results.values()), **results,
            "report_only": {g: rows.get(g) for g in plan["groups"]["report_only"]}}


def grade(plan, sha, runs):
    if set(runs) != set(CONDITIONS):
        raise ValueError("grading requires one result for each of flat, a4 and b4")
    grades, index = {}, None
    for name, run in runs.items():
        design = run["design"]
        case = plan["conditions"][name]
        expected = dict(DESIGN, inputs_sha256=INPUTS_SHA256, genes_simulated=plan["genes"],
                        competitor_transcripts=plan["competitor_transcripts"],
                        prediction_background="reference_index", abundance_ratio=case["abundance_ratio"],
                        seed_offset=case["seed_offset"], replicate_seeds=case["seeds"],
                        g2=_provenance(sha, name))
        for key, value in expected.items():
            if design.get(key) != value:
                raise ValueError(f"{name} result does not match frozen {key}")
        reference = design["reference"]
        if (run["salmon"]["version"] != SALMON_VERSION
                or reference["index"]["salmon_version"] != SALMON_VERSION
                or reference["index"]["k"] != 31
                or not reference["index"]["seq_hash"] or not reference["index"]["name_hash"]
                or not reference["fasta_sha256"]
                or (index is not None and index != reference)):
            raise ValueError(f"{name} result used a different recorded index or Salmon version")
        index = reference
        fragments = run["competitor_fragments"]
        if len(fragments) != DESIGN["replicates"] or any(type(n) is not int or n <= 0 for n in fragments):
            raise ValueError(f"{name} lacks the per-replicate emitted competitor counts")

        rows = {}
        for row in run["genes_tpm"]:
            gene = row["gene"]
            if gene not in plan["genes"] or gene in rows:
                raise ValueError(f"{name} has an unexpected or duplicate TPM row for {gene}")
            n = row["n"]
            if type(n) is not int or not 0 <= n <= DESIGN["replicates"]:
                raise ValueError(f"{name} has an invalid replicate count for {gene}")
            if n >= RULES["minimum_usable_replicates"] and (
                    not math.isfinite(row["bias"]) or not math.isfinite(row["realised_sd"])
                    or row["realised_sd"] < 0):
                raise ValueError(f"{name} has invalid TPM statistics for {gene}")
            rows[gene] = row
        grades[name] = _grade_rows(plan, name, rows)
    return {"passes": all(result["passes"] for result in grades.values()),
            "conditions": grades}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("freeze")
    run = sub.add_parser("run")
    judge = sub.add_parser("grade")
    for p in (prepare, run):
        p.add_argument("--inputs", required=True, type=Path)
        p.add_argument("--background-fasta", required=True, type=Path)
        p.add_argument("--out", required=True, type=Path)
    for p in (run, judge):
        p.add_argument("--plan", required=True, type=Path)
        p.add_argument("--plan-sha256", required=True)
    run.add_argument("--condition", required=True, choices=CONDITIONS)
    run.add_argument("--index", required=True, type=Path)
    run.add_argument("--work", required=True, type=Path)
    run.add_argument("--threads", type=int, default=8)
    for name in CONDITIONS:
        judge.add_argument(f"--{name}", required=True, type=Path)

    a = parser.parse_args(argv)
    try:
        if a.command == "freeze":
            if a.out.exists():
                raise ValueError("choose a new path for the frozen prediction file")
            plan = freeze(a.inputs, a.background_fasta)
            with a.out.open("x", encoding="utf-8") as output:
                json.dump(plan, output, indent=1, sort_keys=True, allow_nan=False)
                output.write("\n")
            print(f"{io.file_sha256(a.out)}  {a.out}")
            return 0
        plan = _read_frozen(a.plan, a.plan_sha256)
        if a.command == "run":
            run_condition(plan, a.plan_sha256, a.condition, a.inputs, a.background_fasta,
                          a.index, a.work, a.out, a.threads)
            return 0
        results = {name: json.loads(getattr(a, name).read_text()) for name in CONDITIONS}
        result = grade(plan, a.plan_sha256, results)
        print(json.dumps(io.finite_json(result), indent=1, allow_nan=False))
        return 0 if result["passes"] else 1
    except (KeyError, OSError, TypeError, ValueError) as error:
        print(f"G2: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
