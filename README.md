# isoform-dominance

[![CI](https://github.com/charliekim97/isoform-dominance-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/charliekim97/isoform-dominance-pipeline/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/isoform-dominance)](https://pypi.org/project/isoform-dominance/)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20672051.svg)](https://doi.org/10.5281/zenodo.20672051)

**Is this functional isoform-class comparison measurable by short reads? Answer that first, then quantify and test it — one config away for any gene.**

Most isoform analyses stumble on two things: deciding *which transcripts form a functional
class*, and knowing whether short reads can resolve the classes at all. `isoform-dominance`
handles both, then quantifies and tests the comparison end-to-end:

```
annotate  gene symbol      → proposed isoform groups        (Ensembl, by 3' terminal exon)
identify  config           → is the class contrast estimable, and how precisely?
extract   Salmon quant.sf  → per-donor isoform-group TPM
stats     per-donor TPM    → paired Wilcoxon, per cohort + stratified combination, + figure
qc        marker TPM       → contamination control (is the signal a cell-type artifact?)
```

![example output](docs/example_output.png)

*Bundled self-test: short LEPR isoform (LepRa) predominates over the long isoform (LepRb) in
control human choroid plexus, in each of two independent cohorts (5/5 and 6/6 donors). The
pooled figure (n = 11, P = 1×10⁻³) is reported alongside the stratified combinations — see
[Statistical notes](#statistical-notes) for why the stratified ones are the numbers to quote.*

---

## Install

```bash
pip install isoform-dominance          # released version
isoform-dominance --version
```

From a clone, for development:

```bash
pip install -e ".[dev]"
pytest -q
```

## Verify it works (no downloads, seconds)

```bash
isoform-dominance selftest
# or: pytest -q
```
reproduces the LEPR reference result (5/5, 6/6, combined n=11 P=9.8e-4) on a clean machine.

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

**2. Identifiability guardrail (`identify`).** A short-read quantifier apportions fragments
by solving a linear inverse problem; a class total is recoverable only when it is an
*estimable* function of that system. This checks the condition directly, on the class-collapsed
system, and reports the smallest fold change the stated design can resolve for each class total
and for the contrast between them. Ask it for the effect size you need:

```bash
isoform-dominance identifiability --config config.json --min-log2fc 0.5
# Identifiability (window=31, k=31, canonical k-mers)
#   annotation: Ensembl release 116
#   background: 9 same-gene transcript(s)
#   design: paired 100bp reads, fragments 200+-60, depth 30M, TPM 10, n=1
#   [identifiable] iso_1165aa: min |log2FC| 0.064; 5369 unique k-mers, 5399 bp in 1 block(s), ~1072 informative reads, conditioning 3.62
#   [identifiable] iso_896aa: min |log2FC| 0.216; 356 unique k-mers, 208 bp in 2 block(s), ~63 informative reads, conditioning 115.05
#   contrast iso_896aa vs iso_1165aa: min |log2FC| 0.234; identifiable, conditioning 117.68
#   effective length: iso_896aa 5103 bp vs iso_1165aa 8072 bp (mean per transcript), log2 ratio -0.66
#   distinguishing windows (0 = 5' end, 1 = 3' end of each transcript): iso_896aa median 0.07 [0.05-0.12], 372 positions over 7 transcript(s); iso_1165aa median 0.67 [0.51-0.84], 10738 positions over 2 transcript(s)
#   unique-read counting floor on log2 ratio: SE 0.187 per donor, min resolvable |log2FC| 0.366 at n=1
#   EFFECT SIZE: |log2FC| 0.500 resolved by both class totals and the contrast at this design
#   VERDICT: identifiable
# exit status 0; with --min-log2fc 0.1 the same run reports NOT resolved and exits 3
```

This example is a counterexample to the skew direction the command states (see below):
iso_896aa is the shorter class but is distinguished at its 5' end (median 0.07), and in the
simulation its bias ran against that direction under 5 of the 6 skews.

The release-116 sequences behind the 2026-09-11 block this replaces, run offline through
`--sequences` and `--background-sequences` with the nine same-gene background transcripts;
transcript and k-mer counts move as the annotation does. The `effective length` and
`distinguishing windows` lines were added from a live run against the same release on
2026-09-13, in which every other line reproduced unchanged. Read the `min |log2FC|` figures first:
they are what the verdict and the exit status are built on. `conditioning` is a diagnostic.
Without `--min-log2fc` the verdict falls back to `--tau` on the conditioning factor, which has
no calibrated value — the same run then reads `weakly_identifiable` — and the command says on
stderr that the verdict is a screening flag, not the exit status.

`effective length` and `distinguishing windows` say how exposed the comparison is to
positional coverage skew, which none of the figures above models (see *It bounds spread, not
accuracy* below). In simulation, Salmon split ambiguous fragments between the classes by
effective length under skew, so the class effective-length ratio tracked the size of the error
and, with the direction of the skew, its sign. Positions are fractions of each transcript's own length, 0 at
the 5' end, because skew acts on each molecule in its own coordinates. When the classes differ
at least 1.23-fold (|log2 ratio| >= 0.3) the command says on stderr which class each direction
of skew tends to inflate, with the evidence: in a 49-gene simulation (Salmon, one quantifier,
monotone positional skew) the 5' direction held for 35–36 of the 39 genes above that ratio, the
3' direction for 30–32, against 20 of 39 at uniform coverage. No other quantifier, gene panel
or form of skew was tested, and below that ratio nothing was measured, so nothing is said.

Exit codes: **0** no `--min-log2fc` given, or it is resolved · **3** `--min-log2fc` given and
not resolved at the stated design · **2** precondition failure: the gene total itself is not
estimable, because a transcript shorter than `--window` has no windows · **1** config or
network error. The structural verdict is never the exit status: it changes with the
annotation release the transcripts came from, so it is reported in the output and in the
`--json` report (`verdict`, with `gene_total` and `effect_resolvable` beside it).

> **Why not just count unique k-mers?** Because that proxy — used by this package up to
> v2.1.1 — is wrong in both directions. A class whose only unique sequence is the ~30 k-mers
> spanning one splice junction *passes* while yielding a handful of usable fragments. A class nested inside another owns no unique k-mer at all and *fails*, yet its
> total is still an estimable function of the system: the containing class is pinned by its own
> unique sequence, so the nested class's indicator lies in the row space. That is an algebraic
> fact about the design matrix — **not** a claim that "the EM will work it out", since an EM
> returns numbers for a non-identifiable model just as readily. Both cases are regression-tested
> in `tests/test_identifiability_model.py`.

**What it checks, in three layers** — all from sequence alone, before any read is quantified:

> **Pass `--background-fasta`.** The default background is the gene's own remaining
> transcripts, which is better than comparing the configured classes alone but still
> blind to pseudogenes, paralogues and homologous loci elsewhere in the index. A
> quantifier resolves fragments against the *whole* index, so any claim about what it
> can separate should be judged against the same FASTA the index was built from. The
> scan is streamed, so a whole-transcriptome background costs memory proportional to
> the gene, not the file.

> **Pin the release.** The verdict is a function of the annotation release: between
> GENCODE v44 and Ensembl 116, 6 of 36 verdicts in a 49-gene panel moved. `annotate`
> records the release it fetched from as `ensembl_release`, and `identifiability` prints it
> in its header and carries it in the `--json` report, beside `fetched_release`, the release
> any sequence was fetched from in that run. `rest.ensembl.org` serves only its current
> release, so a later run that fetches sequence is not a run against the config's release,
> and the command says so when the two differ. `--background-fasta` alone does not pin a
> run: the configured transcripts' cDNA and the gene's other transcripts are still fetched
> live. Pass `--sequences` with the configured transcripts' cDNA and `--background-fasta`
> with the FASTA the index was built from, both from the release the config names, and
> nothing is fetched at all. A verdict quoted without a release is not a reproducible claim.

| Layer | Reports | Why it is not the layer above |
|---|---|---|
| Sequence uniqueness | unique k-mers, bases covered, block structure | judged against a background — the gene's other transcripts by default, the whole index with `--background-fasta` — not just the configured classes |
| Read model | `informative_fraction`, `expected_informative_reads`, counting-noise floor on log2 ratio | informativeness is scored on the **sequenced ends**, not the fragment: a unique region mid-fragment is never observed |
| Estimability | row-space residual, rank, and `min_resolvable_log2fc` from the Poisson-weighted GLS covariance, for each class total **and** their contrast; the structural conditioning factor as a diagnostic | a full-rank system with a near-degenerate contrast passes a rank test and still yields nothing |

The headline figure is `min_resolvable_log2fc`: the smallest |log2 fold change| a 95% interval
excludes zero for at the stated design, under Poisson counting error only.

**It bounds spread, not accuracy.** The figure is built from the delta-method standard error
of the log class ratio, so it describes replicate scatter under uniform coverage and a
correctly specified compatibility model. It says nothing about bias when coverage is not
uniform. In simulation (Salmon, 49 genes, 30 replicates each; 39 genes estimable with a
finite predicted SE in all seven coverage models, each a monotone positional skew), no
gene's |bias| reaches its own predicted SE at uniform coverage: 0 of 39. At a 2.3× ratio
between first- and last-decile gene-body coverage, a routine mildly degraded sample, 20–24
of 39 exceed it, depending on the direction of the skew. The replicate SD meanwhile stays
below the predicted SE in 34–37 of 39, as it did in 35 of 39 at uniform coverage, so the
scatter gives no sign of the bias. NTRK3 holds a replicate SD of 0.016–0.025 while its bias
runs from 0.002 to 2.82 log2. **Replicate agreement does not detect this failure.**
Measure your own libraries' gene-body coverage (RSeQC `geneBody_coverage.py`) before reading
the figure as an error bar.

The conditioning factor is kept as a **diagnostic**, not a headline — a **geometry proxy, not a standard
error** — it assumes `Var(y) = sigma^2 I`, which a quantifier does not satisfy, and it moves by
up to two orders of magnitude with the annotation release. Its thresholds are provisional. The three
layers answer different questions and can disagree. Structural estimability is
not precision, and neither is a promise that a quantifier's optimiser will land on the
right answer — an EM returns numbers for a non-identifiable model too. What the
estimability layer asserts is narrower and checkable: whether the class contrast is a
linear functional of the observable fragment-class expectations.

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

## How this relates to existing tools

Differential transcript usage (DTU) is a mature area, and for genome-wide
discovery you should use the established tools — this one does **not** replace them:

- **DEXSeq, DRIMSeq, satuRn** — genome-wide DTU testing. They assume you already
  have a transcript-by-sample count matrix and defined transcript groups.
- **IsoformSwitchAnalyzeR** — rich functional annotation of isoform switches
  (domains, NMD, coding potential) in R/Bioconductor; grouping and import are
  configured by the analyst.
- **fishpond / swish** — rigorously propagates quantification uncertainty using
  Salmon inferential replicates.

**On identifiability specifically**, this package builds on existing work rather than
claiming the idea:

- **Hiller et al. 2009** ([doi:10.1093/bioinformatics/btp544](https://doi.org/10.1093/bioinformatics/btp544))
  — criteria under which isoform abundances are uniquely determined by RNA-seq observations.
- **Ferrer-Bonsoms et al. 2022** ([doi:10.1093/bioinformatics/btab873](https://doi.org/10.1093/bioinformatics/btab873))
  — a read- and fragment-length-aware identifiability criterion, applied genome wide.
- **terminus** ([doi:10.1093/bioinformatics/btaa448](https://doi.org/10.1093/bioinformatics/btaa448))
  — groups transcripts *post hoc, from the data*, by inferential uncertainty.

The first two address **transcript-level** identifiability of the full deconvolution;
`terminus` chooses its groups after seeing the data. Neither answers the question a biologist
arrives with: *is the contrast between two classes I defined on biological grounds estimable?*

Two different questions hide under the same word, and keeping them apart is the whole reason
this package exists. A contrast is **estimable** when it lies in the row space of the
compatibility system — a yes/no in exact arithmetic. It is **measurable** when it can actually
be recovered at a finite sequencing depth. Transcript-level identifiability settles the first
in one direction only, and is silent on the second.

- **It is not necessary for estimability.** Two transcripts with identical sequence make each
  abundance unidentifiable while their **sum** is perfectly estimable. A transcript-level
  verdict can reject a comparison that is in fact sound.
- **It is sufficient for estimability — and that is the problem.** If every transcript is
  identifiable the design matrix has full column rank, so *every* contrast is estimable. The
  row-space condition is then satisfied trivially and stops discriminating exactly where you
  need it to.
- **It is not sufficient for measurability.** Full rank says nothing about finite depth: the
  class direction can be so weakly observed that a formally estimable contrast is
  unrecoverable. `identifiability` therefore reports the **smallest resolvable effect size
  beside the row-space verdict**, and neither on its own.

`isoform-dominance` targets that question: *for one gene, which functional isoform class
predominates?* The contribution is evaluating the **estimability of the user's own class
contrast, before quantification** — a check none of the tools above performs, and none takes a
bare gene symbol as input. Around that sit two conveniences rather than claims: `annotate` goes
from a gene symbol to a reviewed isoform-group proposal (presented for review, not treated as
final), and the whole thing is a scriptable Python CLI with a download-free self-test, meant to
ship alongside a manuscript.

More documentation: a step-by-step [NTRK2 walkthrough](docs/tutorial_NTRK2.md), a
[gallery of further example genes](docs/example_genes.md) (NTRK2/NTRK3 kinase
truncations and the FLT1 soluble-decoy receptor), and an [API reference](docs/api.md).

## Statistical notes

- Donor-level two-sided **exact Wilcoxon signed-rank** (`scipy.stats.wilcoxon`), per cohort.
- **Small-n floor is computed, not just documented.** `signed_rank_resolution_floor(n)`
  returns `2^(1-n)` — under the sign-permutation null exactly one assignment puts every
  difference on the same side. At n = 5 that is 0.0625, so no arrangement of five donors
  reaches 0.05. Ties among the absolute differences do **not** raise it. Every test is
  reported with its floor and flagged when the floor exceeds 0.05.
- **Cohorts are combined three ways**, all reported: donor-**pooled** (what v2.1 reported alone),
  a weighted **Stouffer** combination of the per-cohort *exact* tests, and a weighted
  **stratified signed-rank** combination. The latter borrows van Elteren's design-free
  `1/(n+1)` stratum weights but is *not* van Elteren's test, which is a stratified
  **two-sample** rank-sum procedure; these data are paired within donor. Stouffer is the
  default, because at these sample sizes each stratum's exact p-value is trustworthy and the
  normal approximation behind a combined rank statistic is not. Pooling donors from independent studies ranks one cohort's differences
  against another's and lets depth or tissue handling drive the result — quote the stratified
  figures when the cohorts are genuinely independent.
- Effect size = median fold-change **with a donor-bootstrap 95% interval** (seeded, so it
  reproduces).

## Layout

```
src/isoform_dominance/   annotate · identifiability · extract · stats · contamination · cli · _selftest
tests/                   pytest (offline; reproduces the reference result + unit tests)
scripts/01_salmon_quant.sbatch
example/                 config + sample maps
.github/workflows/ci.yml docs/  pyproject.toml  CITATION.cff  LICENSE
```

## Citation

**Version provenance.** Version **2.1.1** was used for the LEPR isoform quantification and
aggregation in a leptin-receptor/LRP1 choroid-plexus manuscript (under revision), and is
archived at [10.5281/zenodo.20738150](https://doi.org/10.5281/zenodo.20738150). Releasing a new
version does **not** alter that record: Zenodo mints a separate version DOI and leaves the old
one in place, and the `v2.1.1` tag and release are left untouched by policy — not because they
are technically immutable, but because a published paper cites them. The `extract`
aggregation behaviour those results rest on is unchanged in 2.2.0, and the bundled self-test
still reproduces the same reference numbers.

Cite this repository (see `CITATION.cff`, DOI 10.5281/zenodo.20672051) and Salmon:
Patro, R. et al. *Nat. Methods* **14**, 417–419 (2017). https://doi.org/10.1038/nmeth.4197

## License

MIT (see `LICENSE`).
