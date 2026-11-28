# G2 index-matched abundance simulation

**All three conditions passed the frozen gate** on 2026-11-27. Each of the 30
libraries per condition combined the 49-gene panel with 3,793 additional transcript
origins from the reference FASTA. The additional origins had equal TPM in every
condition.

## Results

For each contrast, bias is the median error of the estimated summed-TPM log2
class ratio relative to the realised fragment and Salmon effective-length truth.
SD is the sample standard deviation of the estimated log2 ratio across usable
replicates. The two metrics are compared separately with the frozen abundance SE
(`min_resolvable_log2fc / 1.96`).

| Condition | Primary \|bias\| < SE | Primary SD < SE | Control \|bias\| < SE | Control SD < SE | Result |
| --- | ---: | ---: | ---: | ---: | --- |
| Required per condition | 6/6 | 5/6 | 37/41 | 32/41 | |
| Flat transcript TPM (`flat`) | 6/6 | 6/6 | 41/41 | 39/41 | Pass |
| Class A:B 4:1 (`a4`) | 6/6 | 6/6 | 41/41 | 39/41 | Pass |
| Class A:B 1:4 (`b4`) | 6/6 | 6/6 | 41/41 | 38/41 | Pass |

No scored contrast reached the stop threshold of absolute bias ≥ 2 × SE.
The control contrasts with SD ≥ SE were CASP9 and COX5A in `flat`, PIK3CA and
PTPRC in `a4`, and ESR1, GHR and PTPRC in `b4`. All primary contrasts had 30
usable replicates. EPOR and PFKL had 25 and 3 usable replicates in `flat` and 28
and 20 in `b4`, respectively. Every other row had 30.

SMN2 and TSC1 had no finite abundance SE in the frozen plan and were report-only.
Both had 30 usable replicates per condition. Values below are in log2 units.

| Gene | Flat bias / SD | 4:1 bias / SD | 1:4 bias / SD |
| --- | ---: | ---: | ---: |
| SMN2 | −0.068 / 0.627 | +2.647 / 0.716 | −0.335 / 0.647 |
| TSC1 | −0.104 / 0.199 | +0.018 / 0.063 | +0.044 / 0.300 |

The pass covers condition-matched predictions using the known simulated
abundances, one depth, uniform coverage and one Salmon index. Calibration for
unknown abundance mixtures, positional bias and real RNA-seq remains open.

## Inputs and outputs

- [Frozen plan](g2-predictions.json) with the 90 seeds, per-contrast predictions,
  eligibility, rules and hashes of the six source files used in the run.
- [Flat TPM](flat.json), [4:1](a4.json) and [1:4](b4.json) run summaries.
  `genes_tpm` holds the scored measurements and `competitor_fragments` counts
  additional-origin fragments in each library.
- [Grade](grade.json) with the per-gene flags and [SHA-256 manifest](SHA256SUMS).

The run used the archived [Ensembl 116 panel](../data/sweep116.json.gz), Python
3.14.6, NumPy 2.5.4 and Salmon 2.7.0 (`-l A`, detected ISF, five threads per
condition). The Ensembl 116 GRCh38 reference-chromosome FASTA contains 446,099
records and has uncompressed SHA-256
`50993c71f1ef9a5b27e2e492bcfba032721c3e784c545a74be6ba02aefa1720a`.
The canonical k=31 index contains 445,269 references after duplicate collapse,
with zero decoys. The summaries record the index's sequence and name hashes.

## Verify the stored grade

From the repository root with the package installed, verify the five files and
recompute the grade using the saved plan. A byte-identical grade makes `cmp`
exit successfully.

```bash
g2_dir=validation/g2
(cd "$g2_dir" && shasum -a 256 -c SHA256SUMS) &&
g2_sha=$(awk '$2 == "g2-predictions.json" {print $1}' "$g2_dir/SHA256SUMS") &&
g2_grade=$(mktemp) &&
python -m validation.gate_g2 grade \
    --plan "$g2_dir/g2-predictions.json" --plan-sha256 "$g2_sha" \
    --flat "$g2_dir/flat.json" --a4 "$g2_dir/a4.json" --b4 "$g2_dir/b4.json" \
    > "$g2_grade" &&
cmp "$g2_grade" "$g2_dir/grade.json"
```

## Rerun the simulation

Use [the Ensembl 116 download](../README.md#remaking-a-design) and the repository's
[reference-chromosome filter](../../README.md#index-scope). Use the resulting
FASTA as `reference_fasta`, with the checkout and NumPy version pinned by the
plan. After the verification commands above, create fresh outputs as follows.

```bash
reference_fasta=Homo_sapiens.GRCh38.cdna.primary.fa
g2_rerun=$(mktemp -d) &&
python validation/run_simulation.py index --fasta "$reference_fasta" \
    --index "$g2_rerun/index" --record-reference --threads 12 &&
(for g2_condition in flat a4 b4; do
    python -m validation.gate_g2 run --condition "$g2_condition" \
        --plan "$g2_dir/g2-predictions.json" --plan-sha256 "$g2_sha" \
        --inputs validation/data/sweep116.json.gz \
        --background-fasta "$reference_fasta" --index "$g2_rerun/index" \
        --work "$g2_rerun/work-$g2_condition" --out "$g2_rerun/$g2_condition.json" \
        --threads 5 || exit
done)
```

Use the same `grade` command with the three new result paths to grade the rerun.
