# isoform-dominance

[![CI](https://github.com/charliekim97/isoform-dominance-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/charliekim97/isoform-dominance-pipeline/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20672051.svg)](https://doi.org/10.5281/zenodo.20672051)

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
#    - from the GENCODE release that matches `annotate --ensembl-release`
#      (https://www.gencodegenes.org/human/releases.html), and
#    - from reference-chromosome transcripts only: see "Index scope" below
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

### Index scope

From GENCODE release 48, `gencode.vX.transcripts.fa.gz` holds the transcripts on
scaffolds, patches and alternate loci as well as those on the reference chromosomes;
releases 44–47 hold the reference chromosomes only. The change is in the release-48 entry
of the changelog in GENCODE's FTP `_README.TXT`; the file description in the same README
still calls the file reference-chromosome only. Ensembl's `cdna.all.fa.gz` holds both too
(checked for release 116). A gene with a copy on such a region (SMN1, HLA-A) is in that
FASTA twice, under its own name and another gene id. An index built from it lets the
quantifier split the gene's reads with the copy — Salmon folds an identical copy into the
reference transcript, but splits reads with a copy that differs, and kallisto splits them
either way — and `extract` counts only the configured transcript ids, so the class totals
are biased. Build the index from reference-chromosome transcripts only. For GENCODE, keep
the transcripts of the reference-chromosome GTF (`gencode.vX.annotation.gtf.gz`):

```bash
gzip -dc gencode.v50.annotation.gtf.gz | awk -F'\t' '$3=="transcript"' \
  | grep -o 'transcript_id "[^"]*"' | cut -d'"' -f2 | sort -u > chr_ids.txt
gzip -dc gencode.v50.transcripts.fa.gz \
  | awk 'NR==FNR {keep[$1]; next} /^>/ {split(substr($0,2), h, "|"); p = (h[1] in keep)} p' chr_ids.txt - \
  > gencode.v50.transcripts.chr.fa
grep -c '^>' gencode.v50.transcripts.chr.fa    # = wc -l < chr_ids.txt  (644,292 for v50)
```

For Ensembl `cdna.all`, keep the records whose header region
(`chromosome:GRCh38:<name>:…`) is 1–22, X, Y or MT:

```bash
gzip -dc Homo_sapiens.GRCh38.cdna.all.fa.gz \
  | awk '/^>/ {split($3, r, ":"); p = (r[3] ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/)} p' \
  > Homo_sapiens.GRCh38.cdna.primary.fa
```

The package checks for such an index where it can. `identifiability
--background-fasta` warns when the FASTA holds a transcript with the configured gene's
name under another gene id, and so does `extract` when `quant.sf` names carry the whole
GENCODE header (an index built without Salmon's `--gencode`). Only an Ensembl header says
where a gene lies, so a same-name gene on a reference chromosome — the chrY copy of a
pseudoautosomal gene such as CD99 or SHOX, or a distinct gene sharing the name — is left
out there, and reported from a GENCODE header as a possible copy. `extract` also writes
`<out>.index.json`, which records the index each donor was quantified against (from
Salmon's `aux_info/meta_info.json`), and it stops when one cohort mixes indexes, unless
`--allow-mixed-index` is given.

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

Cite this repository (see `CITATION.cff`, DOI 10.5281/zenodo.20672051) and Salmon:
Patro, R. et al. *Nat. Methods* **14**, 417–419 (2017). https://doi.org/10.1038/nmeth.4197

## License

MIT (see `LICENSE`).
