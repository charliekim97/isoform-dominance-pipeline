# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses semantic
versioning.

> **Provenance note.** Version 2.1.1 is the version cited for LEPR isoform
> quantification and aggregation in the leptin-receptor/LRP1 choroid-plexus study, and
> is archived at Zenodo DOI 10.5281/zenodo.20738150. That version DOI resolves to its own
> Zenodo record, which a later release does not touch, and the `v2.1.1` tag and release
> are left in place as a matter of policy: they are never retagged, replaced or deleted.
> Nothing in this or any later version alters that record's files. The
> `extract` aggregation behaviour those results depend on — transcript-to-group mapping
> and per-donor TPM summation — is unchanged in 2.2.0, and the bundled self-test still
> reproduces the same reference numbers.

## [Unreleased]

Targets **2.3.0**, not a patch: the exit status of `identifiability` changes meaning.
The version strings in `pyproject.toml`, `__init__.py` and `CITATION.cff` are bumped at
release, once the remaining 2.3.0 items are in.

### Added
- The `identifiability` report carries `gene_total` — the estimability of the sum of every
  column of the compatibility system, with `transcripts_without_windows` listing any
  transcript shorter than `window` — and `effect_resolvable`, which is `true` when both
  class totals and the contrast resolve `--min-log2fc`, `false` when any of them does not
  or has no finite figure, and `null` when no effect size was asked for.
- `identifiability --min-log2fc`: the smallest |log2 fold change| in the class ratio you
  need to resolve. When given it replaces `--tau` in the verdict, and is the recommended
  way to run the command. `--tau` has no calibrated value: across a 49-gene survey the
  median gene sat at conditioning 65.0 against the default tau of 10.0, so the default
  rejects 85% of what it is applied to, and any other fixed value simply sorts genes by
  how many transcripts they have annotated. An effect size is a question the experimenter
  can answer; a conditioning number is not.
- The report now carries, per class and for the contrast, `gls_relative_se` and
  `min_resolvable_log2fc` computed from the Poisson-weighted GLS covariance of the *whole*
  system — the estimator `conditioning_factor` has always described. Under Poisson
  weighting the same 49 genes had a median factor of 16.4 against 65.0 homoskedastic, and
  a 157x range against 850x; the ranking barely moves (Spearman 0.93), so this does not
  rescue a fixed threshold, it puts the number on a scale that converts to an effect.
- `class_coherence` and a `coherence` block per class: the median and minimum pairwise
  window Jaccard among a class's transcripts. `annotate` groups by the 3' terminal-exon
  acceptor alone and asks the user to review the proposal; this is the number to review it
  with. Three of 49 surveyed classes came in under 0.05 — TPI1 at 0.004, whose two members
  are 374 nt and 2217 nt long, CD44 at 0.009, DMD at 0.038 — and the CLI now says so on
  stderr, because a precision figure for such a class describes a quantity nobody asked
  for. Low coherence is a reason to revisit the grouping, not a verdict: most of the
  survey's hard genes had perfectly coherent classes that simply share a lot of sequence.
- `beyond_linear` per estimand, against `LINEARISATION_LIMIT`. Checked against simulation
  (Poisson counts, GLS fit, sample SD of log2(A/B) over 3-4k draws): inside the limit the
  delta-method figure matched at 0.99-1.03x, outside it the two part company by up to
  twofold. Past the limit the figure reads as "not resolvable at this design", not as a
  value.

### Changed
- **Breaking: the structural verdict no longer sets the exit status of
  `identifiability`.** Up to 2.2.0 `not_identifiable` exited 2 and `weakly_identifiable`
  exited 3. The verdict cannot gate a program: it is a function of the annotation
  release — on NTRK3 the same comparison is `not_identifiable` against today's Ensembl and
  `identifiable` against GENCODE v44, and 6 of 36 comparable verdicts in a 49-gene panel
  move between the two — and across seven simulated coverage models the
  `not_identifiable` group never differed from the estimable genes in realised error
  (Mann-Whitney p = 0.08–0.96, with the sign of the difference inconsistent). The exit
  status now answers the question the user asks with `--min-log2fc`:
  **0** when no `--min-log2fc` is given or it is resolved; **3** when it is given and not
  resolved; **2** only when the gene total itself is not estimable, because a transcript
  shorter than `--window` has no windows and an all-zero column — a precondition of the
  system, not a judgement, and independent of the grouping, the design and any threshold;
  **1** for a config or network error, as before. The verdict is still printed and still
  in the `--json` report, and without `--min-log2fc` the command now says on stderr what
  it does and does not mean. A pipeline that branched on exit 2 or 3 to skip
  "unmeasurable" genes should pass `--min-log2fc` with the effect it needs, or read
  `verdict` from `--json`.
- `cli.EXIT_WEAK` is renamed `cli.EXIT_EFFECT_NOT_RESOLVED`; the value is still 3.
- **The conditioning factor is demoted from the headline to a diagnostic** in the CLI
  output, the README and the paper. Each class and contrast line now leads with
  `min |log2FC|`, and the README example is run with `--min-log2fc`. Against realised
  quantification error in simulation, the conditioning factor was the weakest predictor
  measured under both uniform and 5'-skewed coverage (Spearman +0.197 and −0.119, n = 40),
  and it moves by up to two orders of magnitude with the annotation release alone
  (PFKL 2182 → 12.2, TPI1 719 → 7.6, today's Ensembl against GENCODE v44). It is still
  reported, with its `Var(y) = sigma^2 I` caveat, and `--tau` still grades the verdict
  when no `--min-log2fc` is given.
- The CLI no longer prints the counting-noise floor as though it were comparable to the
  per-class figures. It is the precision of `log2(n_a/n_b)` for the informative-read
  counts of each class's *best single transcript*, which equals the class ratio only when
  both classes have the same informative fraction — on TPI1 it reads 0.76 where the class
  totals themselves are not resolvable at all. The line now says which estimand it is.
- Redesigned the `stats` figure. The per-cohort numbers used to float above the axes at
  a fixed offset and collided with the figure title; they now sit in each panel's own
  title block, in reading order (figure title -> combined result -> cohort -> that
  cohort's numbers -> plot). Panels share one y axis, so the slopes are comparable
  across cohorts rather than each panel being independently scaled. Colour now follows
  the isoform **class** and is identical in every panel -- it previously indexed the
  *cohort* on one side of the slope and left the other side grey, so the same class was
  drawn in different hues from panel to panel. The two class colours are a validated
  categorical pair (CVD dE 24.7) and both classes are direct-labelled on the x axis, so
  identity never rests on colour alone. The figure now also shows the bootstrap
  fold-change interval and spells out the exact-test floor when a cohort is
  underpowered.
- The header carries **one** combined P — the default combination, named — instead of
  two competing numbers; the other two combinations are given in a footnote and all
  three remain in the stats CSV. `DEFAULT_COMBINATION` stays Stouffer's, and now says
  why: its inputs are the exact per-cohort tests, so it is never anti-conservative,
  whereas the stratified signed-rank P comes from a normal approximation that is not
  validated at the stratum sizes this package targets. It is reported alongside the
  default, never in place of it.
- Regenerated `docs/example_output.{png,pdf,svg}` and `docs/example_output_stats.csv`,
  which still showed v2.1.1 output (the old 5-column stats table with a single
  `COMBINED` row).
- Added `scripts/make_docs_example.py`, so the committed example figure and stats table
  can be regenerated from the bundled self-test data with no data access. It refuses to
  overwrite the committed files if the reference result does not reproduce. The CSV,
  SVG and PDF it writes are byte-reproducible for a given Matplotlib (`SOURCE_DATE_EPOCH`
  suppresses the embedded timestamp, a fixed `svg.hashsalt` stabilises element ids, and
  the font family is pinned). The PNG is not, across platforms: glyph rasterisation goes
  through FreeType, whose hinting differs by build and CPU architecture, so a
  macOS/arm64 render and a Linux/x86-64 render of the same figure differ in the
  compressed pixel data while looking the same. The committed PNG is the Linux render.
  Previously the committed artefacts had no committed way to rebuild them, which is how
  they came to be a minor version out of date.
- `stats.FONT_STACK` exposes the figure font stack (`["Arial", "DejaVu Sans"]`, unchanged
  by default). Pinning it to one family makes a figure reproduce identically across
  machines; the docs-example script pins DejaVu Sans, which ships with Matplotlib, so
  the committed example is the same file whether it is rebuilt on a machine that has
  Arial or one that does not.
- CI now regenerates the example and fails if `docs/example_output_stats.csv` has
  drifted from what the code produces; a figure difference is reported for review
  rather than failing the build, since fonts and renderers differ across machines.
  `ruff` now covers `scripts/` as well as `src` and `tests`.

### Fixed
- **`scipy>=1.10` was a false floor.** Before 1.15, SciPy's `wilcoxon(method="auto")`
  answered *any* zero difference with the normal approximation, at any n. For the pairs
  in `test_zeros_shrink_the_effective_n_and_raise_the_floor` (five pairs, one tied at
  zero) SciPy 1.14.1 returns P = 0.0455 where the exact value is 0.125: an
  anti-conservative per-cohort p, and the test fails there. The floor is now
  `scipy>=1.15`. The CI matrix installs the newest release of everything, so it could
  never have caught this. A new `floor` job pins every runtime dependency to exactly its
  pyproject minimum and runs the suite and the self-test on Python 3.10. The pins are
  read out of `pyproject.toml`, so the floor and the job cannot drift apart. The
  `wilcoxon_method` code no longer carries the pre-1.15 dispatch rule.
- **Ensembl was asked for one transcript per request, and one failure killed the run.**
  `identifiability` fetched each configured and background transcript with its own GET
  -- 20 to 60 requests for a typical gene -- and none of them was retried. On
  2026-09-11 a degraded Ensembl (HTTP 500/503 and read timeouts) failed every attempt
  at NTRK2, NTRK3 and FLT1. cDNA is now fetched with `POST /sequence/id`, 50 ids per
  request, with results matched back to ids by the echoed `query` field rather than by
  position. Every Ensembl request, including `annotate`'s lookup, is retried on HTTP
  429/500/502/503/504, connection errors and read timeouts. The wait is 1 s before the
  first retry, doubling each time; a 429 waits for its `Retry-After` instead. Other HTTP
  statuses are not retried. Defaults: 5 retries, 1.0 s (up to 31 s of backoff per
  request), 30 s timeout per attempt. The retry count and wait are
  `--retries`/`--retry-wait` on `annotate` and `identifiability`, `retries=`/`retry_wait=`
  in the API, and documented in `isoform_dominance.ensembl`.
- **A background fetch that failed part-way was silently kept.** `analyze` wrapped the
  gene-background fetch in `except Exception` and continued with whatever had arrived.
  A flaky network therefore produced a *different answer* with no error. Observed live:
  a LEPR run reported `background: 1 same-gene transcript(s)` after the lookup had
  succeeded and one cDNA had arrived. Now a gene symbol Ensembl does not know (HTTP
  400/404) still means an empty gene background, and any other failure is raised.
- **A read timeout reached the user as a traceback.** `TimeoutError` is not a
  `URLError`, so it slipped past the `except (URLError, HTTPError)` in `annotate` and
  `identifiability`. It is now reported through the same "Ensembl request failed"
  message, exit 1.
- **README's LEPR `identifiability` example could not be reproduced.** It was captioned
  `--background-fasta gencode.v44… --tpm 5` with nine same-gene background transcripts,
  but its 3399 unique k-mers for `iso_896aa` equal, on the Ensembl release-116 sequences,
  the count with *no* gene background. The silently-partial background fetch above would
  produce exactly that. It is replaced with a default-flag run against release 116,
  2026-09-11: with the complete nine-transcript background, `iso_896aa` has 356 unique
  k-mers, conditioning 115.05, `weakly_identifiable` (exit 3). The drop comes from the
  gene's other transcripts: ENST00001037957 (nonsense-mediated decay) alone takes the
  class to 1125 and `weakly_identifiable`. The README also no longer calls the gallery
  genes "verified short-read separable". Run with defaults on release 116 the same day,
  NTRK2 and NTRK3 are `not_identifiable` and FLT1 is `weakly_identifiable`.
  `docs/example_genes.md`, `docs/tutorial_NTRK2.md` and `docs/index.md` still show the
  earlier output and are not updated here.
- **`selftest` did not accept `--json`**, although the CLI's module docstring (and the
  paper) say every subcommand does; argparse rejected it with `unrecognized arguments`.
  It now writes `{"ok", "checks", "combinations", "headline_combination"}` to stdout.
  Each check names what it compared and the values it saw. Exit codes are unchanged:
  0 pass, 1 fail. `_selftest.result()` returns that object; `_selftest.run()` keeps its
  `(ok, messages)` shape.
- **The self-test reported one of the three cohort combinations.** It printed the
  donor-pooled P only; the Stouffer and stratified signed-rank P existed solely in a CSV
  inside a temporary directory deleted on exit. All three are now printed (and are in
  the JSON), with the headline one marked. They are reported, not checked: the pooled P
  is the reference result and remains the only combination under test.
- **The Wilcoxon p-value's method was hidden.** `wilcoxon` ran with SciPy's default
  `method="auto"`, which is not always exact: a zero difference sends it to the normal
  approximation, and from SciPy 1.15 so does a tie at n > 13 (at n <= 13 a tie selects
  an exhaustive permutation test). Nothing said which. The p-values are unchanged — the
  call still uses `auto`, now named explicitly — and `paired_stat_detail` reports the
  computation that produced them as `wilcoxon_method` (`exact`, `permutation` or
  `asymptotic`). SciPy does not expose that choice, so its rule is mirrored per version
  and the label is kept only if re-running SciPy with the method named explicitly
  reproduces the p-value bit for bit, else `unresolved`. The stats CSV gains a
  `wilcoxon_method` column (appended last, empty on the two combination rows, which are
  not Wilcoxon tests), and the `stats` summary prints it beside each P. Example: 14
  pairs with one zero difference report P = 0.00713 (`asymptotic`); the exact P is
  0.00464.
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
- **`estimability` took its rank and its conditioning factor from two different
  tolerances.** The rank test truncated on `sigma(A)`; the conditioning factor went
  through `pinv(A'A)`, whose `rcond` is relative to `sigma(A)**2`. The two thresholds sit
  five decades apart, and every direction in the band between them was declared estimable
  by the rank test while having its conditioning direction projected away by the
  pseudo-inverse — so the reported factor collapsed to `0.0`, the best possible score, for
  the worst-conditioned contrasts there are. `estimability(diag(1, 1e-7), e2)` returned
  `conditioning_factor = 0.0` where the true value is `1e7`, and `_verdict` graded it
  `identifiable` with no reasons. That inverts the reason the quantity exists. The factor
  is now computed from the SVD the rank test already holds, as
  `sum_i (v_i'c)^2 / sigma_i^2` over the retained directions, so the two cannot disagree;
  a contrast outside the row space reports `inf` rather than a finite number computed on
  its projection. Pinned by `test_conditioning_is_truncated_on_the_same_tolerance_as_the_rank`.
- **`compatibility_matrix` counted distinct window sequences where it needed positions.**
  `A[c, t]` is documented as the probability that a window drawn uniformly from `t` falls
  in class `c`, so each column must sum to 1; with any repeated window it did not. The
  error is one-directional — repetitive shared classes were under-weighted by their repeat
  multiplicity, which is exactly the classes that absorb the most fragments — and repeated
  windows are the rule in cDNA rather than an edge case (A-rich 3' ends, tandem repeats,
  Alu in long UTRs). No estimability verdict moves, because the row space is unaffected,
  but the conditioning factor is a direct function of these numbers. Pinned by
  `test_compatibility_matrix_counts_positions_not_distinct_windows`; the older invariant
  test used a fixture whose windows were all distinct and so could not fail.
- **The zero functional is estimable.** `0'theta = 0` lies in every row space and is
  estimated by the constant `0` with no error; it was reported as `estimable: False`,
  `conditioning_factor: inf`, `rank: 0`. Reachable from a config whose two groups name the
  same transcripts, where the tool refused the self-comparison — correctly — for a reason
  that was not true.
- **Every claim ordering the window system against a real read system is withdrawn.**
  The sentence was wrong four times, in four different forms, each inferring something
  about *row spaces* from something about the *partition*: (1) "anything this system
  declares unidentifiable is unidentifiable for real reads too", stated backwards; (2)
  "anything declared estimable here is estimable from reads", which ignores that reads
  omit; (3) "refining the classes enlarges the row space", which refinement does not
  give; (4) "the per-transcript normalisation re-weights every column, so the rank is
  not monotone", where the normalisation is a positive column scaling and provably
  cannot change the rank at all. What actually holds: a longer window shrinks each
  *position's* compatibility set, but nothing about the resulting matrix is monotone,
  and two distinct things cost the rank rather than one. Four short transcripts now in
  the test suite have surviving-signature counts 4, 7, 4, 3, 5, 5 at windows 3 through 8
  and ranks 4, 4, 3, 3, 4, 4; the contrast is estimable at 3 and 4, not at 5 and 6,
  estimable again at 7 and 8. At window 6 only three signatures survive four
  transcripts, so the count alone caps the rank. At window 5 four survive and the 0/1
  incidence pattern is full rank, yet the matrix is rank three -- the position
  multiplicities are linearly dependent. The normalisation causes neither, and the raw
  count matrix gives the same ranks. The partition of positions also does **not**
  refine: signatures separated at a short window merge at a long one, and the class
  count can fall. The two systems come apart
  both ways: a window no sequenced end can reach contributes a row a real design never
  produces, so a contrast estimable here can be lost; and this construction discards
  adjacency, so two transcripts carrying the same k-mers in a different order are one
  class here while a read spanning the difference separates them. `c in Row(A)` therefore
  neither implies nor is implied by `c in Row(A_observed)` outside `span{1}` -- the
  multiples of the grand total, estimable in any column-stochastic system, `c = 0`
  included. The construction is described as what it is, a sequence-derived screening
  surrogate at a stated window, and whether it predicts what a quantifier recovers is
  handed to the simulation study rather than asserted. Pinned by
  `test_a_longer_window_is_a_different_system_not_a_sharper_one` and
  `test_the_grand_total_is_estimable_in_any_column_stochastic_system`.
- **`--window` is no longer advertised as a sharpness dial.** Its help text said "set to
  the read length for a sharper, still conservative, system" -- both halves of which are
  the withdrawn guarantee, and the first is falsified by the fixture above.
  `_verdict`'s docstring said a functional outside the row space "cannot be recovered at
  any depth"; it now says that of the surrogate system, which is the only thing it is a
  statement about.
- **`informative_fraction` is no longer described as quantifying the omission effect.**
  It reads the group-unique flags and never inspects the shared multi-transcript
  classes, which are also rows of `A`; it quantifies one consequence -- the loss of
  unambiguously assignable fragments -- not the omission of rows.
- **The CLI, README and Statement of Need no longer say "no depth fixes this".** With the
  guarantee withdrawn, exit code 2 cannot claim anything about sequencing depth. It now
  says the contrast lies outside the row space of the compatibility surrogate at this
  window length, and that a longer `--window`, a different grouping or long reads may
  change the verdict.
- **`coverage_stats` double-counted bases where two unique runs were closer than `k`.**
  A run of `r` unique k-mer starts spans `r + k - 1` bases, which is right per run and
  wrong when summed: two runs separated by a gap of fewer than `k` positions have spans
  that overlap or touch. At `k = 3` with unique starts at 0 and 2 the spans are bases
  0-2 and 2-4 — five bases, reported as six. On alternating unique/shared sequence this
  drove `unique_fraction` above 1. `unique_length` and `n_blocks` now come from the union
  of the spans, so a `block` is a maximal run of contiguous unique *bases* — what a read
  has to sit on — rather than a maximal run of unique k-mer starts. `informative_fraction`
  reads the start flags directly and was never affected. Pinned by
  `test_coverage_stats_unions_overlapping_spans`.
- **`sqrt(c'(A'A)^+ c)` was described as the variance factor.** It is its square root, the
  standard-deviation factor; and under `Var(y) = sigma^2 I` there is nothing generalised
  about it, since GLS is OLS there.

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
