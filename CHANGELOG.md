# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses semantic
versioning.

> **Provenance note.** Version 2.1.1 is the version cited for LEPR isoform
> quantification and aggregation in the leptin-receptor/LRP1 choroid-plexus study, and
> is archived at Zenodo DOI 10.5281/zenodo.20738150. That archive is immutable and its
> tag and release are permanent; nothing in this or any later version alters it. The
> `extract` aggregation behaviour those results depend on — transcript-to-group mapping
> and per-donor TPM summation — is unchanged in 2.2.0, and the bundled self-test still
> reproduces the same reference numbers.

## [2.2.0] - 2026-09-09

### Added
- **Group-level estimability.** `identifiability` now builds the fragment-compatibility
  system for the gene and reports whether each class total, and the contrast between the
  two class totals, is an estimable function of it — the textbook row-space condition —
  together with a structural conditioning factor. That factor is `sqrt(c'(A'A)^+c)` and is
  deliberately not called a variance: it is the GLS variance factor under `Var(y)=sigma^2 I`,
  which a quantifier does not satisfy. Its thresholds are provisional. A full-rank system with a
  near-degenerate contrast direction passes a rank test and still yields nothing, so
  both are reported. Formalises at the *class* level what @hiller2009 and
  @ferrerbonsoms2022 established at the transcript level.
- **Background-aware uniqueness.** Uniqueness is judged against the gene's remaining
  transcripts by default (`--no-gene-background` restores the old behaviour), and
  against an arbitrary FASTA via `--background-fasta` — ideally the one the Salmon index
  was built from. The FASTA is streamed, so a whole-transcriptome background costs
  memory proportional to the query rather than the file.
- **Canonical k-mers**, matching what the index actually stores and what an unstranded
  library requires. `--strand-aware` restores the old behaviour.
- **A read/fragment model.** `unique_length`, `unique_fraction`, block structure,
  `informative_fraction` and `expected_informative_reads` at a stated read length,
  fragment-length distribution, depth and class TPM, plus the resulting counting-noise
  floor on the log2 class ratio. Informativeness is evaluated on the *sequenced ends*,
  not the whole fragment: a unique region in the middle of a long fragment is never
  observed.
- **Graded verdicts with reasons**: `identifiable`, `weakly_identifiable`,
  `not_identifiable`, exposed as exit codes 0 / 3 / 2 and as a `reasons` list.
- **`--json` on every subcommand**, emitting the full result object.
- **Statistics**: `min_achievable_p` (the exact-test floor — at n = 5 the two-sided
  floor is 0.0625, so five donors can never reach 0.05), donor-bootstrap confidence
  intervals on the median fold-change, tie accounting and an explicit `zero_method`, and
  two stratified cohort combinations — weighted Stouffer over the per-cohort exact tests
  and van Elteren's design-free stratified signed-rank — reported beside the pooled test.

### Changed
- `identifiability`'s verdict is no longer "does each class own a unique k-mer". That
  proxy erred in both directions and both errors are now regression-tested: a class
  whose only unique sequence is a splice junction passed while being practically
  unmeasurable, and a class nested inside another failed while being perfectly
  estimable. `primary_distinguishable` is retained for callers written against v2.1 but
  is superseded by `verdict`.
- `stats` reports the pooled combination *and* the stratified ones. Pooling donors
  across independent cohorts ranks one cohort's differences against another's; the
  stratified figures are the ones to quote. The `combined` key and the self-test's
  reference numbers are unchanged.
- Package description now leads with the identifiability question rather than with
  quantification.

### Notes
- Public API additions are additive: `paired_stat`, `analyze`, `kmers` and the shape of
  the `stats.run` return value all keep their v2.1 behaviour.

## [Unreleased]

### Fixed
- `identifiability` no longer reports `primary_distinguishable: True` for a
  `primary_comparison` that names a group missing from `groups` (the unknown label
  was silently skipped, so an all-unknown pair reduced to a vacuously true `all()`),
  nor for a single-group config (with no other group, the lone group trivially owns
  every k-mer). Both now raise a clear error. This mattered because the CLI's exit
  code is the gate a pipeline keys on: the bad config exited 0 and `extract`/`stats`
  ran on a comparison that was never checked.
- `identifiability` exit code **1** is new, for a config error. The command now
  distinguishes three outcomes: **0** groups distinguishable, **1** invalid config,
  **2** a primary group has no unique k-mers. Previously a config error surfaced as
  an uncaught traceback, and the two failure kinds were not separable by exit code.

## [2.1.1] - 2026-06-16

### Changed
- Paper and README no longer describe the bundled result as "published"; it is
  framed as a reanalysis of public cohorts (the associated publication citation is
  added on availability), removing an internal inconsistency.

### Fixed
- `qc` validates inputs explicitly: the marker table must contain at least one
  tissue and one contaminant column, the target column must exist, and at least
  three donors must overlap before a Spearman correlation is attempted (previously
  these could silently produce `NaN`).
- `stats` now fails with a clear error when no donors match the requested
  condition, instead of emitting an all-`NaN` figure/table.
- The CLI reports a clear "expected NAME=path" error for malformed
  `--perdonor` / `--markers` / `--target` arguments instead of a raw `ValueError`.
- The Salmon helper script no longer swallows download failures (`|| true`
  removed) and verifies each FASTQ exists before quantification.

## [2.1.0] - 2026-06-14

### Added
- JOSS submission materials: `paper.md` (with Summary, Statement of need, State of
  the field, Software design, Research impact, and AI-usage-disclosure sections)
  and `paper.bib` (incl. NumPy, Matplotlib, Ensembl REST, and dataset citations).
- Community files: `CONTRIBUTING.md` (with a support-channel statement),
  `CODE_OF_CONDUCT.md`, bug-report issue template.
- Documentation: `docs/tutorial_NTRK2.md` walkthrough, `docs/example_genes.md`
  gallery (NTRK2/NTRK3 kinase truncations, FLT1 soluble-decoy receptor — all
  verified short-read separable), and a hand-written `docs/api.md` API reference.
- "How this relates to existing tools" section in the README.
- Unit tests for the CLI (argument parsing and exit codes) and for statistical
  edge cases.

### Changed
- Matplotlib is now imported lazily inside the plotting functions, so
  `--version`, `annotate`, and `identifiability` start without loading the
  plotting stack.
- `paired_stat` handles edge cases explicitly: empty input, all-tied pairs
  (signed-rank test undefined → NaN), and non-finite fold-change ratios.
- `extract` raises clear errors when no `quant.sf` files are found or when a
  file is missing the required `Name`/`TPM` columns.

### Fixed
- Single-group configs (one isoform class, no pair) no longer crash with an
  opaque `IndexError`: `extract` writes the per-donor table without a pair
  fraction, and `primary_pair` raises a clear message when a paired comparison is
  required but only one group is defined.
- `qc` raises a clear error when the config lacks a `contamination_qc` section.
- Removed unused imports; the source now passes `pyflakes` cleanly.

## [2.0.0] - 2026-06-14

### Added
- Packaged, pip-installable `isoform-dominance` with a unified CLI.
- `annotate`: gene symbol → proposed isoform groups from Ensembl (3' terminal-exon clustering).
- `identifiability`: short-read distinguishability check via group-unique k-mers.
- `extract` / `stats` / `qc`: per-donor TPM, paired Wilcoxon, contamination control.
- Offline pytest suite reproducing the published LEPR result; CI on Python 3.10–3.12.
