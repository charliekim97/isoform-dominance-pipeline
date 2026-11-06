# Validation: the simulation behind the documentation's figures (issue #2)

Reads were simulated from 49 genes' real transcript sequence, quantified by Salmon and by kallisto, and the
estimated class ratio compared with the fragments actually drawn, under seven coverage designs. What the package
predicted for each gene, before any read existed, was recorded beside the realised error.

- `data/`: the stored outputs and the inputs they were made from, kept unchanged; `data/README.md` gives the design
  and the SHA-256 of every file.
- `run_simulation.py`, `simulate_reads.py`: the harness that made them. It also records the class ratio of summed TPM
  (below, "Remaking a design"), which the 3.0 gate's run, `data/sim/c1_d1_tpm.json`, holds beside the `NumReads` rows.
- `gate_3_0.py`: grades a TPM run against the gate for 3.0 (below, "The abundance estimand ..."); its grade of the
  stored run is under "The 3.0 gate's result".
- `tpm_run_detail.py`, `whole_index.py`: made the two files computed to describe that run (`data/README.md`, its last
  table): the per-transcript detail of three of its genes, and the 49 genes against the index's FASTA.
- `summarize.py`: recomputes, from `data/` and the installed package, every figure the documentation quotes from
  this simulation, the four-row comparison of issue #2 and the 3.0 gate's result. `summary.json` is its output, and
  the two blocks below are generated from it.

CI does not run Salmon. `tests/test_validation.py` checks that the stored files are the ones `data/README.md` names,
that the package still predicts for every gene what each stored row recorded (when it does not, the stored
comparison no longer describes the code that ships), that `summary.json` and the two blocks below are what
`summarize.py` computes, and that every figure the README, `docs/`, `paper.md`, the CLI, the package's
docstrings, `simulate_reads.py`'s docstring and this file quote from `data/` is the computed one (the comments in
`run_simulation.py` are not checked). `tests/test_gate_run.py` holds the 3.0 gate's run to its grade, the
per-transcript detail to the run's rows, and the whole-index file's gene-background half to what the package gives
now. `tests/test_harness_smoke.py` runs the harness with no quantifier:

- `run` on LEPR and TGFB1, three replicates at `--depth-scale 0.04`, with only the `salmon` command stubbed to write the fragments the simulator drew as `quant.sf`'s `NumReads` (and `meta_info.json`, `lib_format_counts.json`): two rows, each with bias and RMSE 0.
- one replicate of the two genes simulated in two processes, at `PYTHONHASHSEED` 0 and 1: the same truth and, decompressed, the same reads (within one environment; not across numpy versions).

```bash
python validation/summarize.py            # print the summary
python validation/summarize.py --write    # rewrite summary.json and the two blocks below
```

## Remaking a design

Salmon (the stored outputs: Salmon 2.7.0, the bioconda Rust rewrite) and an index of Ensembl 116's cDNA; each
replicate simulates one library of all 49 genes and quantifies it, about 1-2 minutes on a laptop, so a design of 30
replicates takes 30-60 minutes. The harness reads the genes from the stored bundle:

```bash
curl -O https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
python validation/run_simulation.py index --fasta Homo_sapiens.GRCh38.cdna.all.fa.gz --index salmon_idx
python validation/run_simulation.py run --inputs validation/data/sweep116.json.gz --index salmon_idx \
    --work sim --bias-5p 0 --replicates 30 --out c1_d1.json
```

`--bias-5p` is the design (`data/README.md` lists the seven); the seed of replicate r is fixed (`r * 7919 + 13`), so
the same command draws the same fragments. The simulator draws each gene's transcripts in the order the bundle gives
them, which is the order of the per-gene files it was made from (`--rundir` takes those files instead), so the
bundle is read as it is, never sorted; `tests/test_validation_harness.py` checks that the bundle and its unpacked
files give the same replicate. For kallisto, build its index with `index --quantifier kallisto` and run with
`--quantifier kallisto --kallisto-bin <path>`.

Beside the `NumReads` rows (`genes`), a run made by this harness writes `genes_tpm`: per gene, the class ratio of summed
`TPM` against the ratio summed TPM would read had every fragment gone to the transcript it was drawn from -- the
fragments written for each transcript over the effective length the quantifier itself gives it
(`run_simulation.tpm_truth`), a transcript the index dropped as identical to another taking the kept one's -- with the
abundance estimand's prediction. The `NumReads` rows are what they were: `tests/test_harness_smoke.py` runs both with
a stubbed Salmon. `gate_3_0.py` grades a uniform-coverage Salmon run by the gate below, against the predictions fixed
with it, and refuses a run of another design, quantifier or Salmon version, or one whose `NumReads` rows do not
reproduce the stored uniform run (`data/sim/c1_d1.json`, within 0.01 log2), which the same seeds, reads, Salmon and
index give again:

```bash
python validation/run_simulation.py run --inputs validation/data/sweep116.json.gz --index salmon_idx \
    --work sim --bias-5p 0 --replicates 30 --out c1_d1_tpm.json
python validation/gate_3_0.py c1_d1_tpm.json    # exit status 0: passes, 1: fails, 2: not graded
```

## The panel

Every figure below except the row check and the survey, which cover all 49 genes, is over one panel: the genes with a row in all seven Salmon designs, 47 of the 49. A gene has no
row in a design when fewer than three of its 30 replicates gave both classes a positive estimate and a positive
truth; the two left out are listed with the designs they have. kallisto has a row for the same 47 in all seven.
Each figure states its n.

## Results

<!-- summary:start -->
Generated by `python validation/summarize.py --write` from `validation/data`; `tests/test_validation.py` fails when it no longer agrees.

Panel: 47 genes, the ones with a row in all seven Salmon designs (kallisto has a row for the same 47). Left out: EPOR (Salmon rows at bias_5p -1.0, 0, +0.3, +1.0, +2.0); PFKL (Salmon rows at bias_5p 0, +0.3, +2.0).

Every stored row's prediction (verdict, estimable, conditioning, min_resolvable_log2fc, beyond_linear, predicted_se; 720 rows of the 15 files) is what the installed package gives for that gene now: 0 differ.

### The four rows of issue #2

Old gate: the v2.1.1 rule, a strand-aware 31-mer of each class that no other configured class has. Structure: the class contrast estimable at the simulated design. Conditioning: the contrast's conditioning factor against 10, the default `--tau`, which from 2.3.0 grades the verdict only without `--min-log2fc`, and never the exit status. RMSE and excess at uniform coverage (excess = RMSE over the Poisson counting floor of the drawn class sizes); worst |bias|: per gene, the largest of the seven Salmon designs. Medians over the genes of a row.

| issue row | old gate | structure | conditioning | n | RMSE | excess | worst \|bias\| | genes |
|---|---|---|---|---|---|---|---|---|
| 1 | pass | estimable | ≤ 10 | 6 | 0.173 | 1.62 | 2.70 | AR, CASP9, ERBB4, FASN, IGF1R, INSR |
| 2 | pass | estimable | > 10 | 31 | 0.190 | 2.09 | 2.01 |  |
| (none) | pass | not estimable | ∞ | 8 | 0.163 | 6.33 | 2.32 | CANX, GNAS, NTRK2, NTRK3, SMN2, TNFRSF1B, TSC1, VEGFA |
| 3 | fail | estimable | ≤ 10 | 0 |  |  |  |  |
| 3 | fail | estimable | > 10 | 2 | 0.267 | 3.36 | 3.72 | CD44, ESR1 |
| 4 | fail | not estimable | ∞ | 0 |  |  |  |  |

Pass, estimable, conditioning > 10 against ≤ 10 (Mann-Whitney p): uniform RMSE 0.77, uniform excess 0.22, worst |bias| 0.12. ESR1: worst |bias| 6.00; replicates with a defined ratio per design 8/8/19/30/30/30/30. NTRK3's uniform RMSE: 0.017; the estimable genes' median: 0.190.

### S1. Spread, not accuracy (Salmon; the 39 estimable genes with a finite predicted SE)

| design | -2.0 | -1.0 | -0.3 | 0 | +0.3 | +1.0 | +2.0 |
|---|---|---|---|---|---|---|---|
| \|bias\| > predicted SE | 29 | 26 | 20 | 0 | 24 | 33 | 33 |
| replicate SD < predicted SE | 31 | 32 | 34 | 35 | 37 | 36 | 36 |
| median RMSE / predicted SE | 1.94 | 1.78 | 1.18 | 0.61 | 2.54 | 4.14 | 5.03 |

predicted SE = `min_resolvable_log2fc` / 1.96. NTRK3 (not estimable, so not among the 39): replicate SD 0.0165–0.0247, |bias| 0.0016–2.8164 log2 across the seven designs.

### S2. Direction, Salmon (the 39 genes with |log2 effective-length ratio| ≥ 0.3)

| design | -2.0 | -1.0 | -0.3 | 0 | +0.3 | +1.0 | +2.0 |
|---|---|---|---|---|---|---|---|
| 5' rule: sign(bias) = −sign(ratio) | 9 | 7 | 7 | 19 | 35 | 36 | 35 |
| 3' rule: sign(bias) = +sign(ratio) | 30 | 32 | 32 | 20 | 4 | 3 | 4 |

### S7. Direction, kallisto (the same reads) (the 39 genes with |log2 effective-length ratio| ≥ 0.3)

| design | -2.0 | -1.0 | -0.3 | 0 | +0.3 | +1.0 | +2.0 |
|---|---|---|---|---|---|---|---|
| 5' rule: sign(bias) = −sign(ratio) | 10 | 6 | 6 | 7 | 34 | 36 | 35 |
| 3' rule: sign(bias) = +sign(ratio) | 29 | 33 | 33 | 32 | 5 | 3 | 4 |

The 39 include 6 that are not estimable (GNAS, NTRK2, NTRK3, SMN2, TSC1, VEGFA), and are the same genes as those whose ratio exceeds 1.23-fold: yes. kallisto has no bias model and leans to the longer class at uniform coverage too, so its 3' figures are no evidence of the skew; only the 5' ones turn toward it.

### S3. LEPR

log2 effective-length ratio -0.662; its bias ran against the rule in 5 of the 6 skewed designs (bias by design: 0.295, 0.204, 0.087, 0.010, -0.041, -0.055, 0.111).

### S4. Class-length summaries (Salmon, the 47 panel genes)

| summary | Spearman with bias, b = +2.0 | b = −2.0 | b = 0 | 5' sign right at b = +2.0, ratio > 1.23-fold |
|---|---|---|---|---|
| plain mean | -0.547 | +0.169 | -0.180 | 35/39 |
| harmonic mean | -0.541 | +0.164 | -0.181 | 36/41 |
| class total | -0.369 | +0.285 | -0.344 | 35/47 |
| minimum efflen | -0.398 | +0.101 | -0.041 | 32/39 |
| sum of 1/efflen | -0.369 | -0.058 | -0.025 | 24/45 |

"sum of 1/efflen" is log2(S_B/S_A), positive when class A is longer, like the others.

### S5. |bias|, the 8 not estimable against the 39 estimable (Mann-Whitney p)

| design | -2.0 | -1.0 | -0.3 | 0 | +0.3 | +1.0 | +2.0 |
|---|---|---|---|---|---|---|---|
| p | 0.91 | 0.95 | 0.50 | 0.29 | 0.08 | 0.28 | 0.52 |

Not estimable: CANX, GNAS, NTRK2, NTRK3, SMN2, TNFRSF1B, TSC1, VEGFA. Read S8 first.

### S6. What ranks the error (Spearman, the 39 estimable genes)

| design | -2.0 | -1.0 | -0.3 | 0 | +0.3 | +1.0 | +2.0 |
|---|---|---|---|---|---|---|---|
| conditioning vs excess | -0.037 | -0.010 | +0.029 | +0.207 | -0.251 | -0.231 | -0.121 |
| predicted SE vs RMSE | +0.480 | +0.465 | +0.492 | +0.726 | +0.175 | +0.013 | -0.058 |

### S8. The design's confound

Every transcript has TPM 10, so a gene's fragments grow with its transcript count, which is also what makes a class contrast inestimable. Not estimable against estimable, medians: gene fragments 49,980 against 12,966 (3.9x), fragments of the smaller class 2,611 against 298 (8.8x); the counting floor differs at Mann-Whitney p = 0.0013.

### The structural survey of the same 49 genes

Quoted by the package's docstrings, the `--min-log2fc` help and the harness. 1911 transcripts, the shortest 100 nt. The contrast is estimable for 41; their conditioning factor: median 65.0, range 850x, above 10 for 85%. The Poisson-weighted factor at theta = 1 per transcript: median 16.4, range 157x, Spearman 0.93 with the conditioning factor; with the design's expected fragment counts as theta: median 353.74, range 116.8x, Spearman 0.922. Classes whose least similar pair of transcripts shares under 5% of windows: TPI1 iso_63aa 0.004, CD44 iso_536aa 0.009, DMD iso_970aa 0.038; over the 74 classes of the primary comparisons with two or more transcripts, the median share of the least similar pair: 0.49. TPI1's unique-read counting min |log2FC|: 0.76. As the abundance estimand (length-corrected class abundance, `--estimand abundance`) the contrast is estimable for 48: for 7 genes as abundance only (CANX, GNAS, NTRK2, NTRK3, SMN2, TNFRSF1B, VEGFA), and as fragment totals only for none.
<!-- summary:end -->

## What the four-row table shows

Issue #2 asked for four rows: genes the v2.1 unique-k-mer gate passes that the structure calls estimable and well
conditioned (1) or badly conditioned (2), and genes the gate fails that the structure calls estimable (3) or not (4).

- **Row 1** holds 6 genes: AR, CASP9, ERBB4, FASN, IGF1R, INSR.
- **Row 2 did not appear as a separate failure.** The 31 genes with conditioning above 10 had no clearly larger
  error at uniform coverage than the 6 at or below it (median RMSE 0.190 against 0.173, excess 2.09 against 1.62;
  Mann-Whitney p 0.77 and 0.22), and under skew their worst |bias| was smaller (2.01 against 2.70 log2, p 0.12).
- **Row 3** holds 2 genes, CD44 and ESR1, whose median uniform-coverage RMSE, 0.267, is the table's largest (their
  median excess, 3.36, is below the not-estimable row's 6.33). ESR1's worst |bias|, 6.00 log2, is at bias_5p -2.0, where only 8 of its 30 replicates had a defined
  ratio.
- **Row 4** holds no gene: every gene the gate fails is estimable.
- **Outside the four rows**, the 8 genes the gate passes and the structure calls not estimable had no worse RMSE at
  uniform coverage than the estimable genes (median 0.163 against 0.190; NTRK3 0.017). Their large excess (6.33) is
  the confound of S8: as a group they are deeper by construction (a median gene 3.9 times the estimable genes'), so
  most of them have a smaller counting floor; not each one does (SMN2, among them, is the shallowest of all 49 genes).

## Where each figure is quoted

| figures | quoted in |
|---|---|
| S1 | README, "It bounds spread, not accuracy"; `min_resolvable_log2fc`'s docstring |
| S2, S7 | README, the paragraph on `effective length`; the NOTE `identify` prints on the skew direction (`cli.py`); `identifiability.class_efflen_ratio`'s docstring; `DIRECTION_MIN_ABS_LOG2_EFFLEN_RATIO`'s comment |
| S3 | README, after the identifiability example; `class_efflen_ratio`'s docstring |
| S4 | `class_efflen_ratio`'s docstring |
| S5, S6 | CHANGELOG, `[Unreleased]` |
| S8 | `simulate_reads.py`'s docstring; here |
| the table | README, "It bounds spread, not accuracy"; here |
| the survey | the docstrings of `gls_covariance`, `class_coherence`, `min_resolvable_log2fc` and `analyze`; the comment in `_verdict`; the help of `--min-log2fc`; `simulate_reads.py`'s docstring; the README, "Two estimands", and the gate below (the abundance estimand) |
| the 3.0 gate's result | here, the two paragraphs around its block; README, "Two estimands" and "Pass `--background-fasta`"; the NOTE `identifiability` prints on the abundance estimand, its gene-background NOTE when that estimand decides, and the help of `--estimand` (`cli.py`); `identifiability.ESTIMANDS`'s docstring and the comment on `_PRECISION_CHECK`; `docs/api.md`; `docs/tutorial_NTRK2.md`; `paper.md`, "State of the field" and "Software design" |
| D1's result | here, its section; README, "Two estimands", "Pass `--background-fasta`" and "Index scope"; `paper.md`, "State of the field" |

## Read S8 before any comparison of the not-estimable genes with the others

Every transcript was given TPM 10, so a gene's fragment count grows with its number of transcripts, and a large
annotation is also what makes a class contrast inestimable. The not-estimable genes are therefore deeper as a group:
S8 gives a 3.9-fold larger median gene and an 8.8-fold larger smaller class. Any comparison of the two groups (S5,
the table's excess) is confounded by depth. A later simulation should hold the gene's total fixed instead; that
needs the prediction re-solved at the new abundances, which `simulate_reads.py`'s docstring explains.

## The abundance estimand, and what decides whether it becomes the default

Every row here compares a class ratio of `NumReads` with the fragments drawn, and every prediction recorded beside
it is the fragment estimand's: this simulation checks the fragment estimand, which `identifiability` decides by
unless asked. `--estimand abundance` decides by the length-corrected class abundance, whose class ratio is the ratio
that summed TPM estimates, which `extract` and `stats` compare. Whether that estimand is estimable is decided
exactly in the window system `identifiability` builds, and where only it is estimable the answer rests on counting a
transcript's length in window positions (`identifiability.ESTIMANDS`); its precision figures have not been checked
against a quantifier. Whether it becomes the default in 3.0 is decided by the following, fixed on 2026-10-23, before
any such run:

1. **The run.** The simulation is repeated as stored: the same 49 genes, designs, seeds and reads, every
   replicate's library holding all 49 genes, quantified as stored (Salmon 2.7.0 against the stored index,
   `salmon quant -l A` with no option beyond those `run_simulation.py` passes), with the class ratio of summed TPM
   (`quant.sf`'s `TPM`) recorded beside that of `NumReads`. The truth, per replicate, is the class ratio of the sum
   over each class of f_t / E_t, with f_t the fragments the simulator wrote for transcript t (the replicate's
   `fragments_per_transcript`, whose class sums are the truth of the `NumReads` rows) and E_t the `EffectiveLength`
   `quant.sf` gives t in that replicate. A transcript `quant.sf` does not list, because the index kept another of
   identical sequence, takes the E_t of the one kept, and its fragments stay in the truth (within the 49 genes,
   three of VEGFA's `iso_317aa` and three of TP53's `iso_354aa`, each identical to a member of its own class). This
   truth is the ratio summed TPM would read had Salmon given every fragment to the transcript it was drawn from; as
   the fragments drawn do for the `NumReads` rows, it leaves only how Salmon apportions the fragments. The
   simulation's own effective length, max(1, L - 199), is not Salmon's for a transcript shorter than most fragments,
   and with TPM 10 for every transcript as the truth Salmon would be charged for the difference: DMD's 123-nt
   transcript alone would put about 1 log2 between its class's summed TPM and that truth. A replicate counts when
   both classes have a positive estimate and a positive truth. A gene's bias is the median over those replicates of
   the estimated minus the true log2 ratio, its replicate SD the standard deviation of the estimated log2 ratio, as
   in every row here; a gene with fewer than 3 replicates that count fails both criteria and, if it is SMN2 or
   VEGFA, stops the gate as in step 4. For kallisto, reported only, `abundance.tsv`'s `tpm` and `eff_length` stand
   for `quant.sf`'s `TPM` and `EffectiveLength`.
2. **The predictions.** The abundance estimand's contrast `min_resolvable_log2fc` / 1.96 (log2) at the simulated
   design, as this code gives them: `gate_3_0_predictions.json` beside this file, written on 2026-10-23 and kept
   unchanged (SHA-256 `b321a1b10d4319b393801ff0ca53c26e177dcbdff2ae950d76a1976170baae00`). A later change to the
   code does not change them. For the primary set they are CANX 0.337, GNAS 0.117, NTRK2 0.068, NTRK3 0.035, SMN2
   0.442, TNFRSF1B 1.362 and VEGFA 0.501, rounded here; the file holds them unrounded, and the grade is against the
   file. They are figures of the abundance with window positions as lengths, not of summed TPM, which weights a
   transcript by an effective length near L - 199: with L - 199 in place of the window positions none of the 7
   contrasts is estimable in the window system (residual 1e-4 to 5e-2, VEGFA's the largest). The primary set
   therefore tests whether summed TPM meets the figures of a nearby functional.
3. **The primary set, graded first.** The 7 genes of the file's `primary` (CANX, GNAS, NTRK2, NTRK3, SMN2,
   TNFRSF1B, VEGFA: on 2026-10-23 the genes whose contrast is estimable as abundance and not as fragment totals),
   with Salmon, at uniform coverage (`bias_5p` 0). It passes when every one of the 7 has a |bias| below its
   predicted SE, at least 6 of the 7 a replicate SD below it, and step 4 does not stop it: the standard the 39 genes
   of S1 met with `NumReads`, where no |bias| reached the predicted SE and the SD stayed below it for 35.
4. **The stop rule.** Two of the 7 already failed with `NumReads`: SMN2 (|bias| 0.544, SD 0.796, against a
   predicted SE of 0.442) and VEGFA (|bias| 0.549 against 0.501). If either has a |bias| or a replicate SD at or
   above its predicted SE with summed TPM, the default does not change, and why a quantifier's TPM departs from an
   estimable abundance is found first. The one SD step 3 lets a gene miss cannot be theirs.
5. **Predictions past the linearisation limit.** Three of the 7 are past it: SMN2, VEGFA and TNFRSF1B, at a
   relative SE of 0.306, 0.348 and 0.944 on the natural-log scale against the limit of 0.3, where the CLI stars a
   figure and reads it as not resolvable rather than as a value. They are graded by the same first-order SE all the
   same, as S1 graded its 39 genes, 14 of whose predictions are past the limit too. TNFRSF1B's predicted SE, 1.362
   log2, is so wide that its check is a weak one.
6. **The control, graded when the primary set passes.** The 39 genes of the file's `control` (S1's on
   2026-10-23), at uniform coverage, with their abundance predictions from the same file: a |bias| below the
   predicted SE for at least 38, and an SD below it for at least 31, the fewest any design of S1 gave with
   `NumReads`. Against the abundance predictions, `NumReads` at uniform coverage met them for 39 and 34. The control
   checks that the TPM ratio and the abundance prediction agree where the two estimands do.
7. **What decides nothing.** The six skewed designs and kallisto are reported. Under positional skew a gene's
   bias passes its predicted SE often, whichever estimand (S1), and kallisto has no bias model.

If the primary set and the control pass, 3.0 makes `abundance` the default; if the primary set fails, step 4 stops
it or the control fails, 3.0 keeps `fragment`. The run is made once, with the stored seeds. A second run counts only
to correct a fault in the harness found without reading the TPM rows -- a replicate whose `fragments_per_transcript`
differ from those the stored seeds give, a quantification the harness refuses, or a figure not computed as step 1
defines -- with its fix committed before the second run, which is graded in full; both runs are reported. A later
change of the default needs a gate of its own, fixed before its run. `--estimand fragment` stays in 3.0 whatever the
outcome, so that a 2.x report can be reproduced and this simulation checked again.

What a pass can show: in this simulation every transcript has TPM 10, and at that truth the genes whose fragment
contrast is not estimable had no clearly larger `NumReads` |bias| than the others (S5), and 5 of the 7 met the
abundance standard with `NumReads` although their fragment contrast is not estimable. A pass therefore tests the
abundance estimand's precision figures at a flat truth; it does not show that being estimable as abundance is what
makes those 7 comparisons recoverable.

## The 3.0 gate's result

The run was made once, on 2026-10-23 with v2.9.0, as step 1 of the gate says, and is graded below by the gate's own
script: the default does not change. A later change of the default needs, as the gate says, "a gate of its own,
fixed before its run". This gate cannot be graded again by repeating its run: the same seeds give the same reads in
one environment, and the gate's reproduction check (`gate_3_0.check_reproduces`) grades no run whose `NumReads` rows
are not the stored run's, so a second look needs new seeds, a new design or other data. What was computed after the
grade is set apart below it and decides nothing.

<!-- gate-result:start -->
Generated by `python validation/summarize.py --write` from `validation/data`; `tests/test_validation.py` fails when it no longer agrees.

**The grade.** `data/sim/c1_d1_tpm.json` is the run the gate names: its `NumReads` rows reproduce the stored uniform run (largest difference 0 log2). `gate_3_0.py` grades it with exit status 1.

| primary gene | predicted SE | summed TPM \|bias\| | SD | n | \|bias\| below | SD below |
|---|---|---|---|---|---|---|
| CANX | 0.337 | 0.068 | 0.137 | 30 | yes | yes |
| GNAS | 0.117 | 0.002 | 0.069 | 30 | yes | yes |
| NTRK2 | 0.068 | 0.011 | 0.067 | 30 | yes | yes |
| NTRK3 | 0.035 | 0.007 | 0.024 | 30 | yes | yes |
| SMN2 | 0.442 | 0.474 | 0.751 | 29 | no | no |
| TNFRSF1B | 1.362 | 0.072 | 0.314 | 30 | yes | yes |
| VEGFA | 0.501 | 0.024 | 0.269 | 30 | yes | yes |

The |bias| is below the predicted SE for 6 of the 7 and the SD for 6. SMN2's |bias| and SD are at or above its predicted SE, so the stop rule (step 4) stops the change: 3.0 keeps `fragment` as the default. The control is not graded: step 6 grades it only when the primary set passes.

**Beside the grade, deciding nothing.**

- The control, described: a |bias| below the predicted SE for 39 of its 39 genes and an SD below it for 35 (not CSF1R, 0.581 against 0.559; ERBB4, 0.371 against 0.346; INSR, 0.213 against 0.208; PTPRC, 0.126 against 0.078).
- The run's `NumReads` rows against the same predictions: a |bias| below for 5 of the primary set and an SD below for 6; for the control, 39 and 34. VEGFA's `NumReads` |bias|, 0.549, is above its predicted SE, 0.501; with summed TPM it is 0.024.
- LDHA's and NDUFS1's summed-TPM rows (control genes) against a truth that gives the fragments of the class members the index dropped (5 of LDHA's members and 1 of NDUFS1's) to the copies it kept, which no compared class holds (`CHANGELOG.md`, the correction recorded before the run): bias -0.059 and +0.062, against +0.232 and +0.105 with the gate's truth (the run's rows); the SD is the same.
- Three backgrounds (`data/whole_index.json`), computed after the grade. Estimable as abundance and as fragment totals: 48 and 41 of the 49 contrasts with the gene's other transcripts as the background (as every prediction here was made); 47 and 41 against the reference-chromosome records of the index's FASTA, where SMN2's abundance contrast is not estimable, with SMN1's transcripts among its columns (20 SMN1 records on chromosome 5, with 23,037 of the 23,040 windows its competing records share with the configured transcripts); 45 and 40 against the whole FASTA, its records on scaffolds included, where CANX's and PTPRC's abundance contrasts and PTPRC's fragment contrast are not estimable either, their own copies on scaffolds among their columns (19 CANX records on HG30_PATCH, with 20,843 of the 24,123 windows CANX's competing records share with the configured transcripts; 9 PTPRC records on HSCHR1_3_CTG31, with all 9,893 of PTPRC's). None becomes estimable and no predicted SE becomes smaller; the largest increases (log2), the same against both FASTA backgrounds: as abundance TP53 0.194 to 0.241, NR1H3 2.229 to 2.265, ALDOA 0.322 to 0.355; as fragment totals ALDOA 0.309 to 0.354, TP53 0.179 to 0.220, NR1H3 2.394 to 2.431. Beyond those, one verdict changes: CASP9's as fragment totals, `identifiable` to `weakly_identifiable`, against both. Of the 7 contrasts estimable as abundance only with the gene's transcripts as the background, 6 are so against the reference chromosomes (not SMN2's) and 5 against the whole FASTA (not CANX's either).
- Where the reads went (`data/sim/c1_d1_tpm_detail.json`): over the 30 replicates, 245,528 `NumReads`, 0.70% of the 35,032,869 fragments of the libraries, went to records outside the 49 genes; by the gene symbol of those records, LDHA 123,839, NDUFS1 87,700, PTPRC 11,552, CANX 10,472, SMN1 6,967 and SMN2 3,061. 225,428 of them went to 36 records that the index kept in place of transcripts of the 49 genes with the same sequence, so that they hold those transcripts' fragments under another name; 6,829 of SMN1's are among them, on 3 records that stand for 3 of SMN2's transcripts, none of them in a compared class, from which 7,284 fragments were drawn. SMN1's other records were given 138, 0.11% of the 121,779 fragments drawn from SMN2, and SMN2's records outside the 49 genes that stand for none of their transcripts, 2,915.
<!-- gate-result:end -->

**What this shows, and what it does not.** These are findings made after the grade. The predictions the gate fixed
were made, as every prediction of this simulation was, with the gene's other transcripts as the background, which
leaves SMN1 out. The one contrast that failed, SMN2's, is one the package does not call estimable as abundance once
the background is the index's reference-chromosome transcripts. What produced that failure was not measured then: SMN2's
`NumReads` row had already failed before the TPM run (step 4 of the gate). D1, below, measured it: with only the index's
copies of SMN1 and SMN2 on alternate loci left out, SMN1's transcripts on chromosome 5 still in it and no fragment drawn
from SMN1, SMN2's summed TPM met its predicted SE. The finding is also one of scope: a prediction made from the gene
alone can promise what the index cannot give, and the answer to set beside a quantifier's is the one
`--background-fasta` gives with the index's transcript FASTA. Nor does it make the grade a pass: CANX's contrast is not
estimable against the whole index, yet its summed-TPM row was accurate; not estimable means that the contrast is not
determined in general, not that every truth is missed. The simulation's index holds records on scaffolds (alternate
loci, patches), which the README's "Index scope" says to leave out: the copies the index kept in place of LDHA's and
NDUFS1's dropped members are there (`CHANGELOG.md`, the correction recorded before the run), as are CANX's and PTPRC's
copies. The run of step 7, the six skewed designs and kallisto, decides nothing and was not made with this run.

## What produced SMN2's failure (D1), fixed before its run

Step 4 of the gate says that why a quantifier's TPM departs from an estimable abundance is found first. In the gate's
run, SMN2's own records on alternate loci took 2,915 reads, and SMN1's records, besides those that stand for SMN2's
transcripts, 138 ("Where the reads went", above). D1 asks whether the index's records beyond SMN2's own transcripts
produced the failure. It is fixed here, and in `d1_smn2.py`, before any of its libraries is quantified.

- **Reads.** The gate's run's 30 libraries, made again from their seeds with the harness the run used (v2.9.0), and
  checked, replicate by replicate, to hold the fragments the stored truths name.
- **Indexes.** Built as the stored one was (`run_simulation.py index`: Salmon 2.7.0, k 31, no decoys, duplicates not
  kept), from Ensembl 116's `cdna.all` less the records an arm leaves out (`d1_smn2.py fasta`, which refuses another
  input or another result):
  - `ref`: 33 records, the SMN1 and SMN2 records off chromosome 5, their copies on the alternate loci
    HSCHR5_1_CTG1_1 and HSCHR5_2_CTG1_1, which the README's "Index scope" leaves out of an index. The SHA-256 of the
    left-out ids, sorted, one per line, is `e1961e0f15a0f7fd59b1c2a21e8cb1395f7880bb29c52f1cfbf50a0d91e5b434`; of
    the FASTA's text, `697899f6d73e1125808be11bd4ca873f13b185283287c6ac328b82c143d46075`. The index holds 453,527
    records.
  - `gene`: 56 records, every SMN1 record and SMN2's records off chromosome 5, so that of the two genes the index
    holds SMN2's 23 transcripts on chromosome 5 and nothing else, as the prediction did. The SHA-256 of the left-out
    ids is `d62aa9718e9e8644382fe27932a8c1a696c3264bb88eb9642c3c4be522c858fe`; of the FASTA's text,
    `687ae4df6d456a5574951e08b3b68573f329e7761ac429e8125d62f7e3d79572`. The index holds 453,507 records.
- **Quantification.** As the gate's run: each library by Salmon 2.7.0 with the harness's arguments, 12 threads.
- **The grade.** SMN2's summed-TPM row in each arm, against the predicted SE the gate fixed, 0.442
  (`gate_3_0_predictions.json`), as step 3 grades a row: it meets it when its |bias| and its replicate SD are both
  below it. A gene with fewer than 3 replicates that count has no row and does not meet it.
- **What it says.** If `gene` meets it, the failure was the index's: with SMN1 and SMN2 held to SMN2's own
  transcripts, the prediction's scope, SMN2's summed TPM met its figure, and `ref` says whether the README's index
  scope is enough. If `gene` does not, D1 does not put the failure on the index, and what produced it stays
  unmeasured. The `NumReads` rows, and the other 48 genes' rows against the gate's run, are reported and decide
  nothing.
- **What it does not decide.** The default: a later change of it needs a gate of its own, fixed before its run.
- **Once.** Each arm is run once. A second run counts only to correct a fault found before the grade is printed (a
  library without the stored fragments, a `quant.sf` not made against its arm's index, a run that stopped), and both
  runs are reported.

`d1_smn2.py grade` applies this, with exit status 0 when `gene` meets the figure, 1 when it does not, and 2 when the
runs are not the ones named here.

### D1's result

Each arm was run once, on 2026-11-02, and `d1_smn2.py grade` graded `data/sim/d1_ref_tpm.json` and
`data/sim/d1_gene_tpm.json` with exit status 0.

| run | records left out | summed TPM \|bias\| | SD | n | below 0.442 | `NumReads` \|bias\| | SD |
|---|---|---|---|---|---|---|---|
| the gate's run | none | 0.474 | 0.751 | 29 | no | 0.544 | 0.796 |
| `ref` | 33 | 0.110 | 0.352 | 30 | yes | 0.179 | 0.349 |
| `gene` | 56 | 0.112 | 0.353 | 30 | yes | 0.181 | 0.350 |

**The failure was the index's.** With SMN1 and SMN2 held to SMN2's own transcripts, SMN2's summed TPM met its predicted
SE. With only the two genes' copies on alternate loci left out, records the README's "Index scope" also leaves out, it
met it too, although SMN1's 23 transcripts on chromosome 5 stayed in that index, 20 of them among the columns that make
SMN2's contrast not estimable as abundance against the reference chromosomes (above). No fragment is drawn from SMN1 in
this simulation, so those 20 were columns the quantifier could leave empty; the other 3 stand, under SMN1's names, for 3
of SMN2's transcripts of the same sequence, none in a compared class ("Where the reads went", above). D1 does not say
whether the copies would have stopped the gate without SMN1's transcripts on chromosome 5 (no arm leaves out those and
keeps the copies), nor whether SMN2's figure holds when SMN1 is expressed.

**Beside the result, deciding nothing.**

- The other 48 genes' rows, against the gate's run, all moved, by little. Summed TPM: in `ref` the largest bias change
  is 0.0161 (TSC1) and the largest SD change 0.0136 (TNFRSF1B); in `gene`, 0.0161 (TSC1) and 0.1609 (EPOR). `NumReads`:
  the largest bias change is 0.023 (VEGFA) in either arm. EPOR's row counts 24 replicates in `gene` against 25 in the
  gate's run and in `ref`: in one more replicate a class estimate was not positive, so that replicate does not count.
  With the stored index, the gate's run had given the stored uniform run's `NumReads` rows to the last digit; D1 made no
  run with the stored index, so whether these moves come with the indexes or with something else that changed between
  2026-10-23 and 2026-11-02 (macOS 27.0.1 to 27.2; Python, numpy and Salmon the same), it does not say.
- The 3.0 gate's numbers, in either arm: the primary set's |bias| below its predicted SE for 7 of 7 and its SD for 7,
  and the control's for 39 and 35, where the gate asks for all 7, 6, 38 and 31. The gate is not graded again on that:
  it names its run, whose index is the stored one, and a later change of the default needs a gate of its own, fixed
  before its run.

## What these outputs do not cover

- **Depth**: one depth, 30M fragments with every transcript at TPM 10.
- **Class ratio**: the classes as flat TPM gives them; the harness's `--ratio` was not run.
- **Index**: no decoys, Ensembl 116's whole cDNA FASTA, its records on scaffolds (alternate loci, patches) included,
  which the README's "Index scope" leaves out of an index.
- **Genes**: one panel of 49.
- **Quantifiers**: Salmon 2.7.0 and kallisto 0.52.0 only, Salmon with `-l A` and none of its bias corrections
  (`--seqBias`, `--gcBias`, `--posBias`); kallisto has none.
- **Coverage**: monotone positional skew only (`bias_5p`); no GC or fragment-sequence bias, and no other shape.
