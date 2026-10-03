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
> and per-donor TPM summation — is unchanged through 2.5.0, and the bundled self-test still
> reproduces the same reference numbers. 2.4.0 adds two refusals to `extract` where 2.3.0
> wrote a table: a cohort quantified against more than one Salmon index, and a config none
> of whose transcripts is in any donor's `quant.sf`. 2.4.1 counts an index with Salmon's
> decoys and one without as two, and writes the rows in 2.1.1's order again (2.4.0 sorted
> them by donor name).

## [Unreleased]

### Added
- **`annotate --gtf FILE`** (`annotate.build_config_from_gtf`): the groups proposed from a
  local GENCODE GTF -- the comprehensive `gencode.vN.annotation.gtf.gz`, or an Ensembl GTF --
  with no network (issue #13). Ensembl 116 is the last release of the legacy platform, and
  its REST API is kept for 116 only, so an index built from GENCODE 51 or later can be
  matched only by its own files. The config is the one `annotate` writes from REST at the
  release the file is of, `annotation_source` aside: on the 109-gene survey panel, configs
  read from GENCODE 44, 48 and 50 matched REST's at releases 110, 114 and 116 for every gene
  once the reader followed three rules (`scripts/parity/` repeats the comparison; the
  tests hold it on a GENCODE 50 extract of seven genes). A protein's length is
  `(sum of CDS nt + (3 - frame of the 5'-most CDS) % 3) // 3`, which pads a CDS that starts
  mid-codon as Ensembl does (`sum // 3` is wrong for 7,412 transcripts of GENCODE 50); the
  canonical transcript is the one tagged exactly `Ensembl_canonical` (not GENCODE 50's
  `Ensembl_canonical_extended`, and not MANE); and clusters are listed by content (see
  Changed). Ids lose their version, `chr` is dropped and `chrM` is `MT`, and the `_PAR_Y`
  copies of GENCODE 25-43 are left out, as REST before 110 has no chrY gene for them. The
  gene is chosen by the rule REST follows, now `annotate.choose_among`. The release recorded
  is the one the header names (`##description: ... version 50 (Ensembl 116)`);
  `--ensembl-release`, when given, must be the same, and a header that names none is
  accepted with a NOTE. The config gains `annotation_source`: `kind`, `file`, `bytes`,
  `sha256`, `provider`, `gencode_release`, `ensembl_release`, `date`, `description` and
  `n_transcripts` (the gene's transcripts of every biotype). A GENCODE basic GTF -- every
  line below a gene tagged `basic`, in a file of at least 100 transcripts -- is refused: it
  leaves out transcripts the index holds, and at release 110 it changed 51 of 109
  proposals. A symbol no gene has exactly is looked for again ignoring case, as REST does,
  with a NOTE. The file is read in 16 MiB blocks and only a block that holds the gene is
  split into lines.
- **`identifiability --gtf FILE --transcripts-fasta FILE`** (`analyze(..., gtf=,
  transcripts_fasta=)`): the gene's transcripts from the GTF and their sequence from the
  transcript FASTA of the same release, with no network (issue #13). The gene is the config's
  `gene_id`; without one, the one GTF gene that holds every configured transcript; failing
  that, the symbol's, by `annotate`'s rule. It must hold every configured transcript, and
  every transcript of it must be in the FASTA at the GTF's version: a subset FASTA (one of
  GENCODE's basic transcripts lacks 1,917 of the 4,396 transcripts of the 109 survey genes
  at release 116, and `pc_transcripts` every non-coding one) or one of another release
  stops the run with what to pass. The FASTA's ids are read from GENCODE's `|` headers and Ensembl's space-separated
  ones; `_PAR_Y` records are skipped, and an id held twice is an error. The gene
  background is the GTF's every other transcript of the gene, and goes the way a fetched
  one does -- the identical-copy rule, the id merge with `--background-fasta`,
  `--max-window-records`, `--decoys` and `--keep-duplicates` -- so a background the two
  sources give alike gives one report; the tests hold that for LEPR and CD99 against
  recorded REST 116 answers. `--transcripts-fasta` may be the `--background-fasta` too, and
  serves alone a run with `--no-gene-background`. The two replace `--inputs`, `--sequences`
  and `--background-sequences`; `--ensembl-release` must be the GTF's release. The report
  gains `annotation.file_release` and `annotation.source` (the files' names, sizes and
  SHA-256, and what the GTF's header says; None for both from REST), and
  `background.gene_id` is the GTF's gene. `--save-inputs` records the files' release and
  `annotation_source`, and `--inputs` repeats the run. A NOTE says when the config was
  annotated against another release, or proposed from another GTF (by SHA-256).
- **A NOTE when the index is of another release, in both modes**: `--background-fasta`
  holding a configured transcript at another version than the sequence used here, listed in
  `background.fasta_other_versions` (`{id: {"used", "fasta"}}`). The version used is the
  file's, REST's (`ensembl.fetch_cdna_batch(..., versions=)`), or the one a config or
  `--sequences` id carries.
- **With `--gtf`, a same-name gene of the index is placed** (`index_scope.placed_by_gtf`). A
  GENCODE header does not say where a gene lies, so through 2.5.0 a gene of the configured
  gene's name under another id was reported as one that "may be" a copy on a scaffold,
  patch or alternate locus. With the GTF each gains `kind` in `same_name_copies`:
  `reference_gene`, a gene on a reference chromosome -- the chrY copy of a pseudoautosomal
  gene, or HERC3's second gene -- which an index of the reference chromosomes keeps, and
  which is no longer warned of; `off_reference`, placed by a GTF that holds the scaffolds
  and patches; or `not_in_gtf`, which the warning calls a copy off the reference
  chromosomes or a gene of another release.
- **A weekly job holds the GENCODE 50 extract to live REST 116** (`gtf-parity` in
  `ensembl-nightly.yml`): `annotate` for LEPR, CD99 and FOXO1, `identifiability` for LEPR and
  CD99, and whether REST gives the version of the cDNA it serves. An outage is skipped, a
  retired archive and a moved answer warned of; the step is tested against an offline fake.
- **`scripts/parity/`**: `record_rest.py` records REST's answers for a gene list at one
  release, `compare.py` compares them with a GTF and transcript FASTA of that release,
  config by config and cDNA by MD5, and `timing.py` times the file mode. On the extract in
  `tests/data/gencode_mini`: 6/6 configs equal, 217/217 cDNA byte-identical.
- `annotation_files`, which reads a GTF and a transcript FASTA with the standard library
  alone; `io.is_gzip` and `io.open_bytes`.

### Changed
- **`_clusters` and `_proposal.tied_with` are listed by content in both modes**: most
  transcripts first, then the longer representative protein, then the lower terminal-acceptor
  coordinate. Through 2.5.0 clusters with as many transcripts kept the order in which REST
  listed the transcripts, and REST at 116 and the GENCODE 50 GTF list them in different
  orders for 104 of 109 survey genes. A REST config's `_clusters` can be listed in another
  order than 2.5.0 wrote; its groups and `primary_comparison` are unchanged, as the
  proposal already broke every tie by content.

## [2.5.0] - 2026-10-02

### Added
- **`identifiability --decoys FILE`** (`analyze(..., decoys=)`), Salmon's `decoys.txt`, for a
  `--background-fasta` that is a decoy-aware index's gentrome: the records it names are
  skipped unread (issue #3). Every exon is in the genome, so judged against it every window
  inside an exon loses its uniqueness, and the genome record becomes a column of the
  system; but Salmon sets aside only the reads that map better to a decoy than to any
  transcript, so a decoy competes for no read the transcripts explain as well. Without
  `--decoys`, a record longer than 1 Mb (`identifiability.LONG_RECORD`) is listed in
  `background.fasta_long_records` with a NOTE; `--decoys` without `--background-fasta` is
  refused; names the FASTA does not hold are listed in `background.decoys_absent` with a
  NOTE. The CLI and the README now recommend "the transcript FASTA the index was built
  from, without the genome decoys".
- **`identifiability --max-window-records M`** (`analyze(..., max_window_records=)`, default
  20), which bounds the `--background-fasta` records that join the system (issue #14; see
  Fixed).
- **The report names its files by SHA-256**: `background.fasta_sha256` and
  `background.decoys_sha256`, as saved inputs already did for the FASTA. Saved inputs record
  the decoys file too, and a rerun without it, or with another, says so.
- `identifiability.scan_fasta_competitors`, the scan's per-record answer;
  `identifiability.estimability(..., svd=)`; `io.read_decoys`.

### Changed
- **The default background stays the gene's own transcripts, and the NOTE gives the figure
  that decided it** (issue #3). The rule was posted on #3 before the code: had the
  whole-index background changed the exit status at `--min-log2fc 0.5` for 10 or more of the
  100 two-class configs `annotate` proposes at Ensembl 116 (the 109-gene release-series
  panel), `--background-fasta` would have become required. Against GENCODE 50's
  reference-chromosome transcripts (644,292 records) it changed none: every config exits as
  it does without the FASTA. It changed one verdict (SDHA, `weakly_identifiable` to
  `not_identifiable`) and moved the primary contrast's min |log2FC| by a median 0.003
  (q75 0.015; ACO2 the most, 3.8 to 8.6). At the default `--max-window-records 20` a run
  added a median 56.5 records as columns (at most 928) and dropped a median 0.4% of the
  configured transcripts' windows (at most 21%, SDHA; 69 configs lost any), and took a median
  46 s (at most 72 s, four runs at once on four cores) and at most 1.3 GB. Two runs of the
  whole-index arm gave the same answer for every config. In 2.4.1, where the FASTA's records
  were no columns, the same comparison changed no exit status and no verdict.
- **The `--background-fasta` scan is 28 to 45 times faster, with the same answer.**
  Each query window is looked up in both orientations, so a record's windows are no longer
  folded, and a window is read only where a 16-nt seed of it, at a stride of `k - 15`,
  matches the configured transcripts: every matching window holds such a seed, so this is
  exact. A synthetic 210 Mb FASTA took 3.0 s plain and 4.7 s gzipped, against 132 s; a test
  keeps 2.4.1's function and compares the two on 600 random cases. Building and solving the
  compatibility system is faster too (one decomposition for every estimand, no container per
  private window), with figures equal to the last digit; a test keeps the 2.4.1 functions.
- `background` in `--json`: `n_background_transcripts` and `gene_transcripts` count the gene
  background's columns, so a copy of a configured transcript left out is no longer among
  them; `identical_to_configured` is an object whenever `--keep-duplicates` is off (it was
  None without a FASTA); new are `fasta_competitors`, `n_fasta_competitors`,
  `sequence_from_fasta`, `identical_to_background`, `identical_source`, `fasta_sha256`,
  `decoys`, `decoys_sha256`, `decoys_listed`, `decoys_skipped`, `decoys_absent`,
  `fasta_long_records`, `max_window_records`, `fasta_left_out` and `windows_dropped`; and
  `gene_total.transcripts_all_windows_dropped`. The `background:` line counts the FASTA's
  competing records. Saved inputs record `analysis.max_window_records` for a run with a
  FASTA, and a rerun at another M says so.
- **Exit 2 is read from the precondition it names**, a transcript shorter than the window
  (`gene_total.transcripts_without_windows`), not from `gene_total.estimable`. The two agree
  while every column sums to one; with windows dropped by `--max-window-records` a column
  sums to less, and the gene total can be inestimable with no transcript windowless, which
  is no precondition failure. `gene_total.estimable` is still reported.
- Saved inputs record `analysis.keep_duplicates` for a run with a gene background as well as
  for one with a FASTA, since the setting now applies to both; `INPUTS_FORMAT` is unchanged.

### Fixed
The first two entries change an answer.

- **A `--background-fasta` competitor is a column of the compatibility system** (issue #14).
  Through 2.4.1 a FASTA record took windows from the uniqueness layer only, and the system
  was built from the gene background alone; the same competitor passed with
  `--background-sequences` was a column. So the route the CLI recommends gave the
  optimistic answer: on a synthetic two-class case rank 2 and contrast min |log2FC| 0.469,
  exit 0 at `--min-log2fc 0.5`, against rank 3, 0.546 and exit 3. Now a record that holds
  a configured transcript's window found in at most M records of the FASTA
  (`--max-window-records`, default 20) is a column, so a competitor found once gives one
  answer by either route; a record longer than 1 Mb never is. Every window of the system
  that a record left out holds is dropped from every layer -- the unique windows, the read
  model and the system -- for the columns' windows in a first pass over the FASTA and for
  the added records' in a second. The system is then the one with every record a column
  less the rows that touch a left-out record, so its answer is never more optimistic than
  that one's, at any order: a record that shares windows only with a background
  transcript is left out, and the windows it holds go too. A property test compares the
  two on random systems with repeats shared by many records, for M from 0 up. One column
  per sharing record would not run: on a synthetic FASTA in which 3,000 records carried a
  repeat of the configured transcripts it took 2.7 GB, against 90 MB now. A record with
  the id of a gene-background transcript is that transcript, once; when the FASTA holds
  other sequence for it, the column is the FASTA's, with a NOTE.
- **An identical copy in the background is no competitor, wherever it comes from**
  (issue #15). Salmon's default index keeps one of identical sequences. 2.4.1 applied that
  to a FASTA record with a configured transcript's sequence and to nothing else: a
  transcript of the gene background with a configured transcript's sequence took every
  window that transcript had (on a synthetic case, a class with 800 unique windows and
  `identifiable` had none and was `not_identifiable`), and a FASTA record with a
  gene-background transcript's sequence -- the chrY copy of a pseudoautosomal transcript --
  was warned about as a copy that splits reads. The rule
  now covers the whole background: the gene's transcripts fetched from Ensembl,
  `--background-sequences`, the saved background of `--inputs`, and FASTA records. A copy
  of a configured transcript is left out; a copy of another background sequence counts
  once, as the first (gene background by id, then the FASTA in file order); the copy
  warning leaves out both. A NOTE says what a twin means for `extract` when the index keeps
  the twin's name, and how to avoid it. `--keep-duplicates` counts every copy, as before.
- **Only a redirect to Ensembl's page on archives is a retired archive** (issue #12). Any
  redirect that ended outside the REST service was taken for a retirement and not retried,
  so a run that met a maintenance page reported a retired archive. A redirect to
  `www.ensembl.org/help/articles/archives` (or the same path on `ensembl.org`) is
  retirement, with the same message; one that ends anywhere else off the service is
  retried like an outage, and named as one when the retries run out.
- **A `--background-fasta` header whose first `|`-field is empty** (`>|x`, `>|`) is read as a
  record without an id; through 2.4.1 it stopped the run with an `IndexError` traceback.
- **The contrast's `gls_relative_se` and `min_resolvable_log2fc` are infinite when the log
  ratio is not estimable.** The figure is the delta-method SE of log(A/B), a quadratic form
  in its gradient `c_a/a - c_b/b`; through 2.4.1 it was computed whenever the difference
  `c_a - c_b` was estimable, which with class totals that are not leaves the gradient
  outside the row space, and the form is then no standard error (one constructed case read
  0.87 where the system with more rows gave 3.06). The overall verdict and the exit status
  do not change, since both class totals already fail there; the contrast's own verdict
  and reasons can.
- **A compatibility system too large for memory is one line and exit 1**, naming the
  allocation and the likely cause (a genome record, or a repeat many records share), not a
  `MemoryError` traceback.

### Corrections to earlier entries
The entries below are left as written; these sentences in them are wrong.
- 2.4.0, the release series: "8 of 36 verdicts move between GENCODE 44 and release 116 …
  two, CASP9 and NR1H3, gain it. Four of the eight change whether the contrast is estimable
  at all (NR1H3, NTRK3, SMN2, TSC1)". NR1H3's `not_identifiable` at releases 110 to 115 came
  from a gene-background transcript with a configured transcript's sequence, which 2.5.0 no
  longer counts (issue #15); NR1H3 is `weakly_identifiable` at every release. So 7 of 36
  move, six losing resolution and one, CASP9, gaining it, and three of the seven change
  estimability. Of the 62 comparisons fixed at release 116 and posable at all six releases,
  11 change verdict between 110 and 116, not 12.
- 2.4.1, README: "at `--min-log2fc 1.0` it changes across the six for 10 of the 62". After
  2.4.1's own fix of `beyond_linear` it was 12; with 2.5.0 it is 13 (NR1H3 too).

## [2.4.1] - 2026-10-02

### Added
- The author's ORCID in `CITATION.cff`, from which Zenodo takes the creators of each
  release's record.

### Changed
- Files gain fields, and readers of the old ones still work. The `extract` sidecar keeps
  `index_decoy_seq_hash` and `num_decoy_targets` per donor and lists
  `index_decoy_seq_hashes` and `unreadable_meta_info`. Saved inputs record
  `analysis.gene_background` and `sequence_sources`, and `analysis.keep_duplicates` only
  for a run with `--background-fasta`; `INPUTS_FORMAT` is unchanged.
- The `identifiability` header names the window and no longer `k`, and each class line
  reads `class N unique k-mers; best transcript B bp in K block(s), ~R informative reads`.

### Fixed
The first three entries change an answer: each makes the code do what the documentation
already said, or refuses an input that gave a wrong answer. Five more change what 2.4.0
gave for an input it accepted: `extract`'s row order and its decoy rule, an unknown
`gene_id`, a regrouped `--inputs` rerun, and a symbol of another species. The rest stop a
traceback or a misleading message, or fix a record, a script or the documentation.
`extract`'s transcript-to-group mapping and TPM sums are unchanged.

- **`identifiability --min-log2fc`: an estimand past the linearisation limit counts as not
  resolved** (exit 3), as `paper.md` and the 2.3.0 entry below already said it did. Through
  2.4.0 `effect_resolvable` compared the figure with `--min-log2fc` and ignored
  `beyond_linear`, so a contrast at relative SE 0.32, whose first-order figure reads 0.91,
  "resolved" an effect of 1.5 and exited 0. The `EFFECT SIZE` line names the estimands past
  the limit, and the verdict gives it as a reason. The figures and the `beyond_linear` flags
  in `--json` are unchanged.
- **Two configurations that gave a wrong answer are refused** (exit 1, one line).
  A `--window` longer than `--read-length`: no read can hold a whole window, so the
  informative fraction is zero whatever the gene (0.119 at window 100, 0.0 at 101, with
  100-nt reads). A transcript in two groups: in the compared pair its +1 and -1 cancel in
  the contrast, and in any two its column entered the system twice with the first copy all
  zero, so the run reported "the gene total is not estimable", exit 2, naming no
  transcript.
- **`annotate` stops on a human symbol none of whose genes is on a reference chromosome.**
  It proposed groups from an alternate-locus gene and recorded nothing: HLA-DRB3 at release
  116 came from a gene on `CHR_HSCHR6_MHC_APD_CTG1`, and Ensembl 116's `cdna.all` holds 23
  such protein-coding symbols (GSTT1, HLA-DRB4, KIR*, LILRA3, TAS2R45 among them), none of
  which the recommended reference-chromosome index contains. It now exits 1 listing each
  gene of the name with its region (`annotate.NotOnReference`, a `ValueError`); `--gene-id`
  takes one of them anyway, with a NOTE that the recommended index does not contain it.
- **`extract` writes its rows in `quant.sf` path order again**, as 2.3.0 did. 2.4.0 sorted
  by donor name; the two orders part when one donor name is a prefix of another and the
  next character sorts before `/` (`D1`, `D1-2`, `D1.5`), and `stats` bootstraps the fold
  interval by row, so the interval moved (one audited case: [0.365, 1.619] against
  [0.355, 2.204] at seed 7). The transcript-to-group mapping and the TPM sums are
  untouched; a test compares the output byte for byte with a CSV the 2.3.0 code wrote.
- **`extract` counts Salmon's decoys in the mixed-index check.** Two indexes of one
  transcriptome, with and without decoys, share `index_seq_hash` (checked with Salmon
  1.10.3: the decoys go into `index_decoy_seq_hash`, the SHA-256 of nothing when there are
  none). A cohort quantified against both was taken as one index; it is now refused like
  any other mix, and `--allow-mixed-index` combines it with a warning. A donor whose Salmon
  recorded no decoy hash is compared on `index_seq_hash` alone. The sidecar keeps
  `index_decoy_seq_hash` and `num_decoy_targets`, and lists `index_decoy_seq_hashes`.
- **A `meta_info.json` that cannot be read is a warning, not a refusal or a traceback.**
  Empty or cut short it was a one-line refusal the documentation does not list; a JSON list,
  bytes that are not UTF-8 or a hash that is not a string gave a traceback. The donor is now
  treated as having no `meta_info.json`, with one warning naming the file, and recorded
  under `unreadable_meta_info` in the sidecar. `extract` refuses only the two things the
  2.4.0 entry lists.
- **One line and exit 1, not a traceback:** an empty or absent `--quantdir`
  (`extract.NoQuantFiles`, still a `FileNotFoundError`), a sample map with no `donor`
  column, a per-donor table without a class's `_TPM` column, and a `stats --condition` no
  donor has.
- **`extract` names a transcript that two groups share.** Its TPM still goes, as in 2.1.1,
  to the group the config lists last; the warning names the transcript and the groups.
- **`identifiability --inputs` gives a regrouped rerun no gene background its run lacked.**
  The saved file did not say whether the run had one, and a rerun under another grouping put
  every saved transcript the config no longer names into the background. A run saved with
  `--no-gene-background` was then judged against a background the live run never had (one
  audited case: a class's unique count 200 on the rerun, 500 live). Saved inputs record
  `analysis.gene_background`; in a file written before 2.4.1 an empty saved background
  means none. Such a transcript is left out, with a NOTE. `INPUTS_FORMAT` is unchanged:
  2.4.0 reads the new key and ignores it.
- **A `gene_id` Ensembl does not know stops `identifiability`** with one line naming the
  release and the id (exit 1). The HTTP 400/404 fallback meant for a symbol Ensembl does
  not know also swallowed `lookup/id`, so the gene background came out empty and the run
  gave a verdict, exit 0. A symbol Ensembl does not know still leaves it empty.
- **A configured transcript the release lacks is named before the gene background is
  judged.** A config of transcripts only release 116 has, run with `--ensembl-release 110`,
  was refused for a gene background "of another gene"; it now says which cDNA release 110
  does not have.
- **No "reads are split" warning for a same-name copy whose sequence is a configured
  transcript's.** The NOTE beside it says such a `--background-fasta` record is not
  counted, because Salmon's default index keeps one of identical sequences; the warning
  contradicted it, and every pseudoautosomal gene got both against `gencode.v50chr`. The
  warning now leaves out the records the NOTE lists. With `--keep-duplicates` nothing
  changes. `same_name_copies` in `--json` still lists every copy.
- **`identifiability` says which window it used.** Every layer is built from the window,
  and `--k` only sets its default; the header printed `k=15` for `--k 15 --window 31`, and a
  rerun from saved inputs at another k and the same window was told it was a different
  compatibility system. The header names the window, a NOTE says when `--k` went unused,
  and the rerun note compares the window and the k-mer convention only. The per-class line
  separates the class's unique k-mers from its best transcript's stretch, blocks and
  informative reads. `--k` and `--tpm` say in `--help` what they set; `--tpm` is given to
  every transcript of the gene, background included.
- **`identifiability --save-inputs` checks that it can write before the run, and saves
  after the report.** A directory without write permission failed after the run, and the
  failed save came before the report, so the answer was lost with it.
- **Saved inputs say where each sequence came from.** A run with part of its sequence
  supplied (`--sequences`, `--background-sequences`) and the rest fetched stamped the
  fetched release on the whole file. `sequence_sources` gives `"supplied"` or
  `"fetched:<release>"` per id, and `ensembl_release` is recorded only when every sequence
  was fetched from it. A rerun that saves again keeps the sources of the file it ran from.
- **Saved inputs, small things.** A run without `--background-fasta` no longer records
  `keep_duplicates`, so a rerun with a FASTA is not told a record "was not counted then".
  A rerun from inputs that record their release is not called unreproducible for a config
  that records none. A rerun reports the gene id the live run fetched its background as
  (`background.gene_id`), which makes its `--json` report equal the live one apart from
  the `annotation` block, which says where the sequence came from. Saved ids with a version
  (`ENST….3`) are compared without it, as `--sequences` ids are.
- **`annotate` applies the reference-chromosome rule to human only.** The region filter
  and the chrX/chrY rule were applied to every species, so a zebrafish gene on chr23 (or a
  fly gene on 2L, a worm gene on chrI) was dropped for a same-name gene on a chromosome
  with a human name, and nothing was recorded. For another species every gene of the name
  is a candidate, and two or more are `AmbiguousGene` with the `--gene-id` hint.
- **Symbols are quoted into request URLs.** A symbol with a space or a slash made a URL
  `http.client` refuses to send, and that `InvalidURL` was retried as a network failure
  and reported as one. Symbols, species and ids are now quoted as one path segment, and an
  `InvalidURL` is not retried: one line saying the request could not be sent.
- **`annotate --json` keeps its notes**, on stderr: the chrX/chrY choice and the tie that
  chose the alternative class were printed only without `--json`.
- **A JSON file that is not UTF-8 is one line naming it** (config, `--sequences`,
  `--background-sequences`, `--inputs`), not a `UnicodeDecodeError` traceback or, for
  `--inputs`, a "config error" that named no file.
- **A gzipped `--background-fasta` is told by its first two bytes** (`1f 8b`), not by
  its name: a plain file named `.gz` was a `BadGzipFile` traceback, and a gzipped one named
  `.GZ`, or with no extension, a "can't decode byte 0x8b" traceback.
- **`qc` on a config without a complete `contamination_qc` is one line**, naming what is
  missing (`target_group`, `marker_panels`, its `tissue` or `contaminant` list); README step
  4 runs `qc` on the config `annotate` writes, which has none, and that was a traceback. A
  marker or target table without the columns it needs, and a cohort with fewer than three
  donors in both, are one line too (`io.InputError`, a `ValueError`).
- **The weekly Ensembl check tells a retired archive from an outage.** "Ensembl release
  not available", one of the strings it skipped on as a network failure, also begins the
  message for a retired archive, so the release-110 target would have skipped quietly
  every week once that archive is retired. It is now a warning saying the target can no
  longer run. A test runs the workflow's own step scripts against the offline fake.
- **`scripts/01_salmon_quant.sbatch`: one output directory per sample map, and ENA paths
  for every run-number length.** `OUTDIR` defaulted to `$PWD/quant` for every cohort and
  the script skips a donor whose `quant.sf` exists, so the second cohort of the README —
  both example sample maps name their donors `ctrl1`–`ctrl5` — silently kept the first
  one's quantifications. The default is now `$PWD/quant/<sample map name>`, the settings
  can be passed with `sbatch --export`, and the README runs the two cohorts into separate
  directories and extracts both. The ENA directory was built as `0` and the last two
  digits, right only for 8-digit run numbers (`ERR2060213` is under `ERR206/003`, not
  `013`); it now follows ENA's rule for 6 to 9 digits and stops on anything else. The
  README gives the `salmon index` command for the filtered FASTA.
- **Documentation that said what the code does not do.** The Wilcoxon test is SciPy's
  `method="auto"`, exact only without zeros or ties (README, `docs/api.md`, `stats.py`, which
  called it exact). Only the per-cohort and pooled tests carry a floor (README and `stats.py`
  said every test). A `--background-fasta` scan holds one record at a time, so its memory is
  set by the longest record (one 20 Mb record: about 110 MB resident), not by "the gene" or
  "the query". `docs/api.md`: zero unique k-mers does not mean a class cannot be measured;
  `kmers` folds to canonical by default; `stats.run` returns six keys, not two; the import
  line lists `ensembl`. README: DEXSeq tests exon usage from its own exon-bin counts;
  Kmerator also takes gene names, so that is not a claim; the v2.1.1 record is cited by the
  manuscript under revision, not by a published paper. CONTRIBUTING lists the CI jobs.

### Corrections to earlier entries
The entries below are left as written; these sentences in them are wrong.
- 2.4.0, saved inputs: "the rerun's whole `--json` report now equals the live one" — it
  differed in `background.gene_id` (fixed above) and differs, by design, in the
  `annotation` block, which says where the sequence came from.
- 2.3.0, `DEFAULT_COMBINATION`: Stouffer's "inputs are the exact per-cohort tests, so it is
  never anti-conservative". A per-cohort p is exact only without zero differences or ties;
  otherwise it is a permutation p (n ≤ 13) or the normal approximation, as
  `wilcoxon_method` reports, and the guarantee does not hold.
- 2.3.0, the README LEPR example: "`weakly_identifiable` (exit 3)". With the default flags
  2.3.0 exits 0 there: from 2.3.0 the verdict does not set the exit status, and exit 3
  needs `--min-log2fc`.
- 2.3.0, "`identifiability` exit code **1** is new … **0** groups distinguishable, **1**
  invalid config, **2** a primary group has no unique k-mers", and "exit code 2 cannot
  claim anything about sequencing depth": these describe the exit status before 2.3.0. From
  2.3.0 it is 0 / 3 (`--min-log2fc` not resolved) / 2 (the gene total not estimable) / 1,
  as that entry's "Breaking" item says.
- 2.2.0, statistics: `min_achievable_p` — the key is `p_floor`, and the column
  `resolution_floor_P`.

## [2.4.0] - 2026-09-30

### Added
- **`--ensembl-release N` on `annotate` and `identifiability`** (`release=` in
  `annotate.build_config`/`run`, `ensembl_release=` in `identifiability.analyze`): fetch
  the annotation, the configured transcripts' cDNA and the gene background from Ensembl
  release N instead of the current one. This closes issue #5, which 2.3.0 only half-closed
  by recording the release: a config's release can now be fetched again, and when a run's
  release differs from the config's the note names the flag that fetches the config's.
  `rest.ensembl.org` serves only its current release. Earlier ones are served by Ensembl's
  REST archive through an alias per release, `e110.rest.ensembl.org`, which answers with
  an HTTP 301 to a date-named host, `jul2023.rest.ensembl.org`. `ensembl.resolve_server`
  follows that redirect once, with a GET, and every later request goes to the host it
  lands on, because `urllib` replays a POST that receives a 301, 302 or 303 as a GET with
  no body: against the alias, the batched cDNA fetch fails with
  `400 {"error":"ID '' not found"}` (observed 2026-09-24). The host must report exactly N on `/info/data` before anything is
  fetched from it. The alias is asked first, so a pinned run against an archived release
  never touches `rest.ensembl.org`, the server whose outage on 2026-09-11 is recorded
  below; only when the alias fails is `rest.ensembl.org` asked which release is current.
  The current release has no working alias (`e116` answered 503, and later redirected to
  a host that answers 503) and is taken from `rest.ensembl.org`. A release later than the
  current one, an archive that reports another release, an archive Ensembl has retired
  and an archive still unreachable once the retries are spent are one-line errors and
  exit 1 (`ensembl.ReleaseNotServed`, a `LookupError`). A retired archive is named as
  such at once, without retries, when its alias redirects off the REST service to a web
  page about archives, as releases 90 to 104 but 94 did on 2026-09-24; releases 75 to 89
  and 94 answered 503 instead, and an archive that does not answer is either down or
  retired, which the message says. Release 111's archive timed out on most lookups
  through 2026-09-24, and 106's intermittently. A 200 that is
  not a release listing, such as a maintenance page, counts as an unreachable archive
  rather than a config error. `--ensembl-release` takes a positive integer. When the
  pinned release has no cDNA for a transcript the config names, the error names the
  release and the transcript; a batch `POST /sequence/id` returns the ids it knows and
  silently drops the rest, so the missing transcript surfaces only there. A run on
  supplied sequence resolves nothing and makes no request, as before. Every request takes
  `server=`, the base URL, which defaults to `ensembl.SERVER`.
- **`identifiability --save-inputs FILE` and `--inputs FILE`**: a verdict that outlives
  Ensembl's REST service. `--save-inputs` writes the sequence the run used — the configured
  transcripts' cDNA, the gene background after the configured transcripts are removed from
  it, and the release they came from — to a JSON file (`"format":
  "isoform-dominance/inputs/1"`); `--inputs` repeats the run from that file with no request
  at all. Ensembl 116 (June 2026) is the last release on the REST API, archives are
  retired as they age (none from 104 or earlier answered on 2026-09-24), and one can be
  down for hours (111 that day), so the file, not the archive, is what keeps a verdict
  reproducible. `--inputs` replaces `--sequences`, `--background-sequences` and
  `--ensembl-release`, and refuses a file saved for another gene, or one that lacks a
  transcript the config names rather than fetching it: a silent top-up would mix releases.
  A `--background-fasta` is not copied into the file, but its path, size and SHA-256 are,
  and a rerun given no FASTA, or another one, says so on stderr. A `--save-inputs` path in
  a directory that does not exist, or that is itself a directory, is refused before
  anything is fetched, not after. The
  header names the saved release, and the `--json` report carries it as
  `annotation.inputs_release`. **The file holds the gene, not one grouping of it.** Which
  transcripts are a class and which are background is a property of the config, not of the
  file, so a rerun pools the saved `sequences` and `background_sequences` and splits them
  again by the config in hand — the rule a live run follows. Splitting them by the saved
  grouping instead lost a transcript the config no longer named from the background as
  well as from its class, and refused one moved the other way for "no sequence" while
  holding that sequence. On LEPR at release 116 with `--min-log2fc 0.3`, against the
  release-116 proposal with one transcript dropped from `iso_896aa`: live gives a contrast
  of 0.323 and exit 3, and the saved-grouping rerun gave 0.264 and exit 0 — the optimistic
  answer, with no warning. The rerun's whole `--json` report now equals the live one,
  for that config, for the full proposal and for the reviewed two-class config of
  `docs/api.md`. The file also records the `k`, the window and the k-mer convention the
  run used, none of which is in the config, and a rerun at another one says so on stderr
  as it does for a changed `--background-fasta`. Every field `--inputs` goes on to read is
  checked when the file is loaded, so a hand-edited or truncated file is one line naming
  the file rather than a `KeyError` from inside the analysis. In the
  library, `analyze(..., inputs_out={})` hands back what was used, and
  `io.save_inputs`/`io.load_inputs` write and read the file.
- **`annotate` flags a proposal decided by a tie.** The alternative class is the
  non-canonical cluster with the most transcripts, and a tie goes to the longer
  representative protein — a tie-break, not a biological criterion. When another cluster
  has as many transcripts as the one chosen, `annotate` now says so on stderr, naming each
  tied cluster by protein length and terminal acceptor, and the config records the rule
  and the tied clusters under `_proposal` (`alternative_rule`, `tied_with`). The note
  quotes `alternative_rule` as recorded rather than paraphrasing its first step: a tie on
  transcript count is settled by protein length and then by the terminal-acceptor
  coordinate, and a note that stops at the protein does not say what decided this one. Across six
  releases of a 109-gene survey, 26 to 39 of the 84 to 100 two-class proposals were such
  ties. Of the 84 genes proposed a pair at every release, 30 were proposed a different
  pair at some release, always on the alternative side; 12 of the 30 had been a tie at the
  release before the first change, and in 8, LEPR and NTRK3 among them, the new
  alternative was a cluster the old one had been tied with. 33 of the 34 changes came
  with releases 115 and 116. Re-proposing all 654 stored configs of that survey from
  their stored lookups with this version gives the same groups, pairs and clusters for
  652; FOXO1 at 115 and 116 changes, through the tie-break fix below.
- **`extract` records which index quantified each donor, and refuses a cohort that mixes
  indexes (issue #8).** Salmon writes `<donor>/aux_info/meta_info.json` with the index's
  `index_seq_hash`, `index_name_hash`, `num_valid_targets`, `keep_duplicates` and
  `salmon_version`. `extract` now copies these to a sidecar `<out>.index.json`, not to the
  per-donor CSV, which `stats` reads, and stops before writing anything when the donors
  of one cohort have different `index_seq_hash` values
  (`index_scope.MixedIndexError`, exit 1): until now it summed TPMs from different indexes
  without a word. `--allow-mixed-index` (`allow_mixed_index=True`) combines them anyway,
  with a warning, and the sidecar says `"mixed": true`. A donor with no `meta_info.json`
  (kallisto, or Salmon output without `aux_info`) is a warning.
- **A same-name copy of the gene on a non-reference region is reported (issue #8).** From
  GENCODE release 48, `transcripts.fa.gz` holds the transcripts on scaffolds, patches and
  alternate loci as well as those on the reference chromosomes; releases 44–47 hold the
  reference chromosomes only. The change is recorded in the release-48 entry of the
  changelog in GENCODE's FTP `_README.TXT`, whose file description still calls the file
  reference-chromosome only. Ensembl's `cdna.all` holds both too (checked for release 116).
  An index built from such a FASTA lets the quantifier split a gene's reads with its
  copies, which `extract` does not count. `identifiability --background-fasta` now lists,
  under `background.same_name_copies`, every gene id in the FASTA that has the configured
  transcripts' gene name but not their gene id, and warns on stderr; `extract` does the
  same when `quant.sf` names are whole GENCODE headers (an index built without Salmon's
  `--gencode`). Only an Ensembl header says where a gene lies, so only there is a same-name
  gene on a reference chromosome (1–22, X, Y, MT) left out; from a GENCODE header it is
  reported as a possible copy, and the warning says what else it may be: in GENCODE v50's
  reference-chromosome FASTA, 34 names with a protein-coding transcript have more than one
  gene id, among them the chrY copies of pseudoautosomal genes such as CD99 and SHOX, which
  GENCODE 44, 48 and 50 give their own gene id, and distinct genes sharing a name such as
  HERC3. The copies are read from the FASTA headers (GENCODE fields 2 and 6; Ensembl
  `gene:`, `gene_symbol:` and the region), not from Ensembl REST:
  `xrefs/symbol/homo_sapiens/SMN1` also returns SMN2's gene ids. On GENCODE v50
  `transcripts.fa.gz` this reports one copy for SMN1 (`ENSG00000275349`) and seven for
  HLA-A, and none for LEPR, NTRK2, NTRK3 or FLT1, whose documented numbers are therefore
  unaffected. The README's "Full workflow" step 0 and `scripts/01_salmon_quant.sbatch` now
  say to build the index from the GENCODE release matching `--ensembl-release` and from
  reference-chromosome transcripts only, with a filter for each source.

### Changed
- A release after 116 is refused as not on the REST API, not as a release that "does not
  exist yet": Ensembl 116 is the final release of the legacy platform, whose REST API
  "remains available for e116 for long term use" with "no plans to port over" to the new
  platform (ensembl.info, 2026-07-21). Unpinned runs therefore read release 116 for as
  long as the service lasts (`ensembl.REST_LAST_RELEASE`).
- The weekly Ensembl check (`ensembl-nightly.yml`) runs a second target, LEPR pinned to
  release 110 through the REST archive, against its own recorded answer
  (`.github/lepr-e110-baseline.json`). An archived release is frozen, so a move there is a
  change in this package or in the archive, and the archive path has its own ways to
  break. Workflows use `actions/checkout@v5` and `actions/setup-python@v6`, which run on
  Node 24.
- Until this release the version on `main` was `2.4.0.dev0`, so that a saved-inputs file
  and a `--version` from an unreleased tree could not be mistaken for 2.3.0's.

### Fixed
- **`identifiability --background-fasta` counted a record with a configured transcript's
  sequence as competing sequence (issue #10).** Such a record took every unique k-mer the
  transcript had, and Salmon's index does not hold it: unless built with
  `--keepDuplicates`, it keeps only the first of identical sequences. In GENCODE v50's
  reference-chromosome index Salmon removed 1,600 (644,292 in the FASTA, 642,692 targets),
  among them all 382 chrY transcripts of the 18 protein-coding genes on both chrX and chrY;
  214 of the 1,600 were identical to a transcript of their own gene. A CD99 config of the
  chrX gene (42 transcripts, each byte-identical to a chrY record) had no unique k-mer in
  either class against a FASTA of CD99's chrX and chrY records, and 1,273 and 3,051 without
  it. Against the whole reference-chromosome `transcripts.fa`, with the gene's other
  transcripts as background too, it had none either (`weakly_identifiable`); it now has 615
  and 2,546 and is `identifiable`, and LEPR's example config gives the same answer both
  ways. A background FASTA record whose whole sequence, upper-cased and stripped, is a
  configured transcript's is now left out, and reported under
  `background.identical_to_configured` (record id -> the configured transcript) and in one
  line on stderr. `--keep-duplicates` (`keep_duplicates=True`) counts it as before, for an
  index built with `--keepDuplicates`; `--save-inputs` records the setting under
  `analysis`, and a rerun from the file with the other setting is noted. Two configured
  transcripts with one sequence are left as they were. The gene background fetched from
  Ensembl is not changed: a transcript there identical to a configured one of the same gene
  is still counted.
- **`annotate` chose silently between genes that share a symbol (issue #9).** Ensembl's
  `lookup/symbol` returns one gene per symbol and does not say when there are others: for
  the pseudoautosomal genes CD99, CRLF2, CSF2RA, IL3RA and SHOX it returns the chrY copy
  (releases 110 and 116), and for HERC3 (110, 116) and DUSP13B (116) the newer of two genes
  of that name on one chromosome, with 1 and 2 transcripts against 22 and 8. All 42
  transcripts of `annotate --gene CD99` belonged to ENSG00000292348, on chrY.
  `build_config` now also asks `xrefs/symbol` for the genes of the name and keeps, with one
  `lookup/id`, those whose display name is the symbol and that lie on 1–22, X, Y or MT. Of
  a chrX/chrY pair it takes the chrX gene and says so; on any other set it stops
  (`annotate.AmbiguousGene`, exit 1) and lists them. `--gene-id` (`gene_id=`) names the
  gene outright. The config records `gene_id`, and `_gene_choice` when there was a choice.
  The chrX gene is the copy a GENCODE-built Salmon index keeps: Salmon keeps only the first
  of identical sequences, GENCODE's FASTA lists chrX first and Ensembl's cDNA FASTA chrY
  first. On a 132-record index of CD99's records, Salmon 2.7.0 kept chrX in GENCODE v50
  order and chrY in Ensembl 116 order. A symbol with one gene costs one request more than
  before (`xrefs/symbol`), the only request that lists every gene of a name.
- **`extract` said nothing when `quant.sf` lacked the transcripts the config names
  (issue #9).** In the GENCODE v50 reference-chromosome index Salmon removed 1,600
  transcripts as duplicates, among them all 66 chrY transcripts of CD99, and the CD99
  config above came out as 0.0000 TPM and a ratio of NA for every donor, with exit 0. A
  cohort in which no donor's `quant.sf` has any configured transcript is now refused before
  anything is written (`index_scope.NoConfiguredTranscripts`, exit 1). Some missing is one
  warning for the cohort, with the count, the first five ids, the donors and the two usual
  causes, and `missing_transcripts` in `<out>.index.json`.
- **`identifiability` fetched the gene background by symbol (issue #9).** For a config of
  CD99's chrX gene, that background was the chrY gene's 66 transcripts, whose sequence is
  the configured one: no unique k-mer in either class and `not_identifiable`, against
  `identifiable` (1,273 and 3,051 unique k-mers) with the same sequence and no gene
  background. A config with `gene_id` is now looked up by id; one without it is looked up
  by symbol as before, and refused when that gene holds none of the configured transcripts.
  The report gives `background.gene_id`, `--save-inputs` records `gene_id`, and `--inputs`
  refuses a config whose `gene_id` is another gene's.
- **An archive alias that redirects to `rest.ensembl.org` itself was read as a retired
  archive.** The test for "still inside the REST service" accepted only hosts ending in
  `.rest.ensembl.org`, and `rest.ensembl.org` is not a subdomain of itself, so
  `--ensembl-release 116` would have failed with "the REST archive for Ensembl release 116
  is retired" the day Ensembl pointed `e116.rest.ensembl.org` at the current server. 116 is
  the last release on the REST API and so the release everyone will pin, and its alias has
  no archive of its own: on 2026-09-24 it answered 503, then redirected to a host that
  answers 503. `ensembl.resolve_server` now accepts the host of `ensembl.SERVER` as well as
  its subdomains.
- `annotate` printed a traceback for a gene with no protein-coding transcript; it now
  prints one line and exits 1.
- **Network failures raised while a response was being read escaped the retry loop.**
  `urllib` wraps an error raised while *sending* a request in `URLError`, but not one
  raised while reading the answer, so these arrived raw — measured against a local socket
  server: a connection closed before the status line (`RemoteDisconnected`), a garbled
  status line (`BadStatusLine`), a body shorter than its `Content-Length`
  (`IncompleteRead`), a reset mid-body (`ConnectionResetError`). They were neither retried
  nor reported: one attempt, then a traceback. They are now retried like any other
  transient failure (`ensembl.TRANSIENT`), and once the retries are spent `annotate` and
  `identifiability` print one line and exit 1. A 200 whose body is not JSON, such as a
  maintenance page, is retried too (`ensembl.BadResponse`), for the batched cDNA fetch as
  well as for lookups (`ensembl.request_json`).
- **`annotate`'s proposal depended on the order in which the server listed transcripts**
  when two alternative clusters tied on both transcript count and protein length: the sort
  was stable, so the first one listed won. Nothing guarantees that order. 249 of the
  14,054 two-class genes of GENCODE 50 (1.8%) have such a tie, and for FOXO1 and STK11
  the REST order and the GTF's order pick different alternatives. Ties now go to the
  lower terminal-acceptor coordinate, and a gene without a canonical transcript takes the
  longest protein, then the lower acceptor, so the proposal is a function of the
  annotation alone; `ALTERNATIVE_RULE` names the new last step. Of the release series'
  654 stored configs, FOXO1 at 115 and 116 changes; no number quoted here or in the
  README moves (the release report and its independent recount are unchanged).
- A file named on the command line that is missing, unreadable or not what it should be
  — a config that is not JSON, not an object or has no `groups`; a `--sequences` or
  `--background-sequences` file that is not an object of id to sequence — printed a
  traceback. Each is now one line naming the file, and exit 1 (`io.InputError`,
  `io.load_json`). `qc` still accepts a config with no `groups`.
- **The annotation-release numbers in the 2.3.0 entry came from a reconstruction of
  GENCODE 44, not from GENCODE 44.** It restricted release 116's transcripts to the ids
  GENCODE 44 lists, which kept release 116's sequences. Rebuilt from release 110 as served
  (release 110 is GENCODE 44), through `--ensembl-release`, on the same 36 genes: the
  transcript sets were right to two transcripts in 797, but 154 of the 795 it kept have a
  different sequence at release 110, and two of 36 verdicts were wrong. NR1H3 is
  `not_identifiable` at release 110, not `weakly_identifiable`: there a class transcript
  and a background transcript have identical cDNA, which puts a null direction across the
  class boundary, and release 116 lengthened one of them by 300 nt. CASP9 is
  `weakly_identifiable` at release 110, not `identifiable`: two of its class transcripts
  were 2.6 and 3.1 kb shorter, which leaves one class 36 expected informative reads
  against the floor of 50 and the contrast a conditioning factor of 11.6 against `--tau`
  10. So **8 of 36 verdicts move between GENCODE 44 and release 116, not 6**, and not all
  in one direction: six lose resolution and two, CASP9 and NR1H3, gain it. Four of the
  eight change whether the contrast is estimable at all (NR1H3, NTRK3, SMN2, TSC1); the
  other four only cross `--tau`, which has no calibrated value, and two of them the
  informative-read floor as well. The median fold change of the contrast's conditioning
  factor is 3.2 on the same genes; 2.3.0's analysis quoted 4.6, which cannot be recovered
  from the two result files it was computed from: that median is 3.8 there, and none of
  some 770 other summaries of them (per class, pooled, with and without the non-estimable
  genes, medians of ratios and ratios of medians) gives 4.6. What stands: 960 of 1911
  transcripts across 49 genes are new since GENCODE 44, 47 of 49 genes changed, NTRK3's
  comparison is `identifiable` at GENCODE 44 and `not_identifiable` at release 116, PFKL's
  conditioning is 12.15 against 2182 and TPI1's 7.62 against 719. Six releases, 110 and
  112 to 116, put a date on it: the transcript count of a 109-gene survey goes 2005, 2062,
  2062, 2066, 3188, 4396, so a comparison measured at any release from 110 to 114 was
  measured against nearly the same annotation, and 115 and 116 more than doubled it. The
  growth is GENCODE's long-read annotation: every one of the 1122 transcripts release 115
  added to those genes, and 1128 of the 1211 that 116 added, come from TAGENE, GENCODE's
  manually supervised pipeline for long-read transcriptome data (Mudge et al., *Nucleic
  Acids Res.* 2025, doi:10.1093/nar/gkae1078). Release 115's were all protein-coding;
  639 of 116's are nonsense-mediated-decay candidates. Genome-wide, GENCODE lists 89,843
  protein-coding transcripts at release 48 (Ensembl 114), 211,446 at 49 (115) and 278,455
  at 50 (116), and 21,902, 21,949 and 91,818 nonsense-mediated-decay ones; Ensembl put
  the release-115 addition at about 121,000 protein-coding transcripts. Release 111's
  archive did not answer long enough to be read.

## [2.3.0] — 2026-09-15

A minor version, not a patch: the exit status of `identifiability` changes meaning.

### Added
- The `identifiability` report carries `gene_total` — the estimability of the sum of every
  column of the compatibility system, with `transcripts_without_windows` listing any
  transcript shorter than `window` — and `effect_resolvable`, which is `true` when both
  class totals and the contrast resolve `--min-log2fc`, `false` when any of them does not
  or has no finite figure, and `null` when no effect size was asked for.
- `identifiability --min-log2fc`: the smallest |log2 fold change| you need both class totals
  and the contrast to resolve. When given it replaces `--tau` in the verdict, and is the recommended
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
- **The annotation release is recorded and reported (issue #5).** `annotate` writes
  `ensembl_release`, taken from Ensembl's `/info/data`. `identifiability` prints it in its
  header and carries it in the `--json` report as `annotation.ensembl_release`, beside
  `annotation.fetched_release`: the release any sequence was fetched from in that run, and
  `null` when none was. The two are different facts. `rest.ensembl.org` serves only its
  current release, so a config annotated against one release and re-run later can fetch
  sequence from another, and the command then says so on stderr. A config without a
  release still runs; the command prints one line saying the verdict is not reproducible.
  The verdict is a function of the release: between GENCODE v44 and Ensembl 116, 960 of
  1911 transcripts across 49 genes are new, 47 of the 49 genes changed, and 6 of 36
  verdicts moved. `--background-fasta` alone does not pin a run, because the configured
  transcripts' cDNA and the gene background are still fetched live; `--sequences` with
  `--background-fasta`, both from the config's release, makes no request at all.
- **The report says how exposed a comparison is to coverage skew.** The contrast carries
  `log2_efflen_ratio`, the log2 ratio of the classes' plain mean effective lengths
  (`max(1, L - frag_mean + 1)` per transcript), with `class_mean_efflen`, and
  `distinguishing_window_position`: median and quartiles of where each class's
  distinguishing windows start, as a fraction of each transcript's own length. The CLI prints
  both. When `|log2_efflen_ratio| >= 0.3` (`efflen_direction_in_band`) it also says on stderr
  which class a 5'- or 3'-skewed library tends to inflate, with the evidence: in a 49-gene
  simulation (Salmon, one quantifier, monotone positional skew) the 5' direction held for
  35–36 of the 39 genes with a ratio above 1.23x and the 3' direction for 30–32, against 20
  of 39 at uniform coverage. Below the band nothing was measured, so nothing is said. Five
  summaries of class length were scored on that simulation: plain mean (Spearman −0.547 at
  b = +2.0, +0.169 at b = −2.0, sign 35/39), harmonic mean (−0.541, +0.164, 36/41), class
  total (−0.369, +0.285, 35/47), minimum (−0.398, +0.101, 32/39), sum of 1/efflen (−0.369,
  −0.058, 24/45). The plain mean is best or tied. The class total weights by transcript count,
  reads −0.344 at b = 0 where the answer has to be null, and is excluded for that reason. The
  table is kept in the `class_efflen_ratio` docstring.

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
- `annotate` no longer writes `"reference": "Ensembl REST (live annotation)"`. The string
  named no release, and `ensembl_release` replaces it. Nothing read the field, so configs
  that carry it load unchanged.
- **`min_resolvable_log2fc` is documented as a bound on spread, not on accuracy.** It is
  built from the delta-method SE of the log class ratio under the Poisson-GLS covariance,
  which assumes uniform coverage and a correctly specified compatibility model. In
  simulation (49 genes, 30 replicates, 39 estimable with a finite predicted SE in all seven
  coverage models), |bias| exceeds the predicted SE in 0 of 39 genes at uniform coverage
  and in 20–24 of 39 at a 2.3× first/last-decile coverage ratio, while the replicate SD
  stays below the predicted SE in 34–37 of 39 (35 of 39 at uniform). NTRK3 holds an SD of
  0.016–0.025 while its bias runs from 0.002 to 2.82. The README says, next to the figure,
  that replicate agreement does not detect this failure.

### Fixed
- The `identifiability` step of the weekly `Ensembl live check` workflow
  (`ensembl-nightly.yml`) would have failed at its first scheduled run: it runs LEPR against
  live Ensembl with default flags, where LEPR is `weakly_identifiable` (release 116), and
  that exited 3. The workflow was added in 2.2.0 and had not yet run. The exit-status change
  above resolves it.
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
