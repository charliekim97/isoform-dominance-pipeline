# isoform-dominance

[![CI](https://github.com/charliekim97/isoform-dominance-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/charliekim97/isoform-dominance-pipeline/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20692725.svg)](https://doi.org/10.5281/zenodo.20692725)

**Isoform-usage quantification *and discrimination* from bulk RNA-seq — one config away for any gene.**

Most isoform analyses stumble on two things: deciding *which transcripts form a functional
group*, and knowing whether short reads can even *tell the groups apart*. `isoform-dominance`
handles both, then quantifies and tests the comparison end-to-end:

```
annotate  gene symbol      → proposed isoform groups        (Ensembl, by 3' terminal exon)
identify  config           → are the groups distinguishable by short reads? (k-mer uniqueness)
extract   Salmon quant.sf  → per-donor isoform-group TPM
stats     per-donor TPM    → paired Wilcoxon, per cohort + combined, + figure
qc        marker TPM       → contamination control (is the signal a cell-type artifact?)
```

![example output](docs/example_output.png)

*Bundled self-test: short LEPR isoform (LepRa) predominates over the long isoform (LepRb) in
control human choroid plexus across two independent cohorts; combined n = 11, P = 1×10⁻³.*

---

## Install

```bash
pip install -e ".[dev]"        # from a clone
isoform-dominance --version
```

## Verify it works (no downloads, seconds)

```bash
isoform-dominance selftest
# or: pytest -q
```
reproduces the published LEPR result (5/5, 6/6, combined n=11 P=9.8e-4) on a clean machine.

## What makes it more than a quantifier

**1. Auto isoform grouping (`annotate`).** Give a gene symbol; it pulls the gene's
protein-coding transcripts from Ensembl, clusters them by their 3' terminal-exon splice
acceptor (the alternative last exon that defines functional isoform classes), and proposes a
comparison you review and rename.
```bash
isoform-dominance annotate --gene LEPR --out config.json
#   iso_896aa : 7 transcripts   (short / LepRa)
#   iso_1165aa: 2 transcripts   (long / LepRb, canonical)
```

**2. Identifiability guardrail (`identify`).** Short reads can only quantify an isoform group
that has *unique* sequence. This checks per-group unique k-mers and **refuses to pretend** a
group is measurable when it isn't.
```bash
isoform-dominance identifiability --config config.json
#   [OK] iso_896aa : 3399 unique k-mers
#   [OK] iso_1165aa: 5369 unique k-mers
#   primary_comparison distinguishable by short reads: True
```

## Full workflow (real data)

```bash
# 0) build a decoy-aware index once (Salmon + GENCODE) — see scripts/01_salmon_quant.sbatch
# 1) quantify on an HPC cluster:
sbatch scripts/01_salmon_quant.sbatch                    # -> quant/<donor>/quant.sf
# 2) extract per cohort:
isoform-dominance extract --config config.json --quantdir quant \
    --samplemap example/sample_map_GSE228458.csv --cohort GSE228458 --out perdonor_GSE228458.csv
# 3) stats + figure:
isoform-dominance stats --config config.json --condition control \
    --perdonor GSE228458=perdonor_GSE228458.csv --perdonor GSE137619=perdonor_GSE137619.csv \
    --out results/dominance
# 4) optional contamination control:
isoform-dominance qc --config config.json \
    --markers GSE228458=markers_228.csv --target GSE228458=perdonor_GSE228458.csv --out results/qc
```

## Statistical notes

- Donor-level two-sided **exact Wilcoxon signed-rank** (`scipy.stats.wilcoxon`), per cohort + combined.
- **Small-n floor:** n = 5 cannot reach P < 0.05 in the exact two-sided test (floor 0.0625);
  report exact P + direction (e.g. 5/5) and combine concordant cohorts for the summary statistic.
- Effect size = median fold-change, reported alongside significance.

## Layout

```
src/isoform_dominance/   annotate · identifiability · extract · stats · contamination · cli · _selftest
tests/                   pytest (offline; reproduces the published result + unit tests)
scripts/01_salmon_quant.sbatch
example/                 config + sample maps
.github/workflows/ci.yml docs/  pyproject.toml  CITATION.cff  LICENSE
```

## Citation

Cite this repository (see `CITATION.cff`, DOI 10.5281/zenodo.20692725) and Salmon:
Patro, R. et al. *Nat. Methods* **14**, 417–419 (2017). https://doi.org/10.1038/nmeth.4197

## License

MIT (see `LICENSE`).
