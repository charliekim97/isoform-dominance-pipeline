# API reference

`isoform-dominance` is primarily a command-line tool, but every subcommand is a
thin wrapper over a small public Python API. Each module is importable from the
`isoform_dominance` package and can be used directly in scripts or notebooks.

```python
from isoform_dominance import (annotate, contamination, ensembl, extract, identifiability,
                               index_scope, io, stats)
```

The config object passed throughout is a plain dict (loaded from JSON via
`io.load_config`) with at least:

```json
{
  "gene": "LEPR",
  "groups": {"short": ["ENST00000371060", "ENST00000616738"], "long": ["ENST00000349533"]},
  "primary_comparison": ["short", "long"]
}
```

`annotate` also writes `ensembl_release`, the Ensembl release the groups were proposed
from. It is optional, so configs written before 2.3.0 still load, but `identifiability`
reports a verdict without it as not reproducible.

---

## `annotate`

Propose isoform groups for a gene from Ensembl by clustering protein-coding
transcripts on their 3' terminal-exon splice acceptor.

**`annotate.fetch_transcripts(gene, species="homo_sapiens", server=None) -> dict`**
Query Ensembl for the gene `lookup/symbol` gives and return `{gene, gene_id, species, strand, transcripts: [{id, protein_aa, terminal_acceptor, is_canonical}]}` (`annotate.transcripts_of(record, gene)` does the same for an expanded lookup record). Raises `ValueError` if no protein-coding transcripts with a translation are found. `server` is a base URL from `ensembl.resolve_server`; None is the current release. Network access required.

**`annotate.cluster_by_terminal_exon(info) -> list`**
Group the transcripts from `fetch_transcripts` by terminal-exon acceptor coordinate; returns clusters, each `{acceptor, rep_aa, n, canonical, ids}`, listed by content: most transcripts first, then the longer representative protein, then the lower acceptor coordinate. So are a config's `_clusters` and `_proposal.tied_with`; through 2.5.0 clusters with as many transcripts were listed in the order the source gave the transcripts, which REST and a GTF give differently.

**`annotate.propose_groups(info) -> (groups, primary, clusters)`**
Choose the canonical cluster and the largest alternative cluster, returning `groups` ({label: [ids]}, labels like `iso_896aa`), the `primary_comparison` list (alternative first), and all clusters. The rule is `annotate.ALTERNATIVE_RULE`: the non-canonical cluster with the most transcripts, a tie going to the longer representative protein, and a tie on both to the lower terminal-acceptor coordinate, so the proposal is a function of the annotation alone.

**`annotate.alternative_ties(clusters, groups, primary) -> list`**
The clusters the alternative class was chosen over on the tie-breaks rather than on transcript count: those, other than the two proposed, with as many transcripts as the alternative. Empty when the choice was not a tie. A release that adds one transcript to a tied cluster changes the proposal; across six releases of a 109-gene survey, 26 to 39 of 84 to 100 proposals were ties.

**`annotate.choose_gene(gene, species, looked_up, **net) -> (record, choice)`**
The gene a symbol means. `looked_up` is the expanded record `lookup/symbol` gave; `xrefs/symbol` lists the other genes of the name and one `lookup/id` keeps those whose display name is the symbol and that lie on 1–22, X, Y or MT (`xrefs/symbol/SMN1` also lists SMN2 and alternate-locus copies). One such gene: `(it, None)`, with no `lookup/id` when xrefs lists no other gene. A chrX/chrY pair — a pseudoautosomal gene — gives the chrX gene, the copy a GENCODE-built Salmon index keeps (`annotate.GENE_RULE`). Any other two or more raise `annotate.AmbiguousGene` (a `ValueError`) listing each gene's id, location and transcript count. `choice` is `{rule, chosen, candidates, reason}` when there was a choice to make, or when xrefs listed no gene and the others were not looked for.

**`annotate.choose_among(gene, named, species="homo_sapiens") -> (record | None, choice | None)`**
The rule both sources follow, `choose_gene` from REST and `build_config_from_gtf` from a GTF: among `named`, the expanded records of every gene whose display name is the symbol, the human candidates are those on a reference chromosome; one is taken with no choice recorded, of a chrX/chrY pair the chrX gene, and any other two or more raise `AmbiguousGene`. `NotOnReference` when every one of `named` lies off the reference chromosomes.

`annotate.NotOnReference` (a `ValueError`) is raised for a human symbol none of whose genes lies on a reference chromosome (1-22, X, Y, MT); the message lists each gene with its region, and `gene_id=` takes one of them anyway. The region and chrX/chrY rules apply to `homo_sapiens` only (`annotate.REFERENCE_SPECIES`); for any other species two genes of the name are `AmbiguousGene`.

**`annotate.build_config(gene, species="homo_sapiens", release=None, gene_id=None) -> dict`**
Convenience wrapper returning a complete, reviewable config dict (including `gene_id`, the Ensembl gene the groups were proposed from, `ensembl_release`, `_proposed` notes, `_clusters`, and `_proposal`: `alternative_rule` and `tied_with`, the clusters from `alternative_ties`, and `_gene_choice` from `choose_gene` when there was a choice). `gene_id` names the gene outright, by `lookup/id`, and must be a gene of the symbol; it is `--gene-id`. `release` is the Ensembl release to propose the groups from; None is the one `rest.ensembl.org` currently serves, and an earlier one is read from Ensembl's REST archive (see `ensembl.resolve_server`). `ensembl_release` records the release the server reported.

**`annotate.build_config_from_gtf(gene, gtf, species="homo_sapiens", gene_id=None, release=None, notes=None, block=annotation_files.BLOCK) -> dict`**
`build_config` from a local GTF -- GENCODE's comprehensive `gencode.vN.annotation.gtf.gz`, or an Ensembl GTF -- with no network (issue #13). The GTF is read by `annotation_files.scan` as REST-shaped records and the gene chosen by `choose_among`, so the config is the one `build_config` writes from the release the file is of, but for `annotation_source`: `{kind: "gtf", file, bytes, sha256, provider, gencode_release, ensembl_release, date, description, n_transcripts}`, `n_transcripts` the gene's transcripts of every biotype. `ensembl_release` is the release the header's `##description` names (`... version 50 (Ensembl 116)`); `release`, when given, must be the same, or `ValueError`, and is recorded as the file's when the header names none. `gene_id` must be a gene of the GTF named `gene`. A symbol no gene has exactly is looked for again ignoring case. A GENCODE basic GTF raises `annotation_files.BasicGTF`. `notes`, a list, receives what a REST run cannot see: a header with no release, a symbol found only ignoring case. `block` is the size of the pieces the file is read in.

**`annotate.run(gene, out, species="homo_sapiens", release=None, gene_id=None, gtf=None, notes=None) -> dict`**
As `build_config`, or `build_config_from_gtf` when `gtf` is given, but also writes the config JSON to `out`. Backs the `annotate` CLI subcommand; `release` is `--ensembl-release`, `gene_id` is `--gene-id`, `gtf` is `--gtf`.

```python
cfg = annotate.build_config("FLT1")
cfg["groups"]               # {'iso_1338aa': [...], 'iso_733aa': [...]}
cfg["primary_comparison"]   # ['iso_733aa', 'iso_1338aa']
```

---

## `identifiability`

Check whether isoform groups can be distinguished by short reads.

**`identifiability.fetch_cdna(transcript_id, **retry) -> str`**
Fetch one transcript's cDNA sequence from Ensembl. Network access required.

**`ensembl.fetch_cdna_batch(ids, versions=None, **retry) -> dict`**
`{id: cdna}` via `POST /sequence/id`, 50 ids per request, matched back to ids by the
`query` field Ensembl echoes. `analyze` fetches all of its sequence this way. `versions`, a
dict, receives each sequence's versioned id when the answer gives one.

**`ensembl.fetch_release(**retry) -> int`**
The release the server is serving, from `GET /info/data`; raises `ValueError` unless it
lists exactly one. `annotate.build_config` records it, and `analyze` reports it as
`annotation.fetched_release` when, and only when, the run fetched sequence.

**`ensembl.resolve_server(release=None, **retry) -> str`**
The REST base URL that serves Ensembl `release`, for the `server=` argument of every other
function here. None, or the current release, is `ensembl.SERVER`
(`https://rest.ensembl.org`). An earlier release is served by Ensembl's REST archive
through an alias, `https://e<N>.rest.ensembl.org`, that redirects (HTTP 301) to a
date-named host -- release 110 to `https://jul2023.rest.ensembl.org`. The redirect is
followed once, with a GET, and the host it lands on is returned; the alias itself is never
used for the batched cDNA fetch, because `urllib` replays a redirected POST as a GET with
no body and Ensembl answers that with `400 {"error":"ID '' not found"}`. The host must
report exactly `release` on `/info/data`. The alias is asked first, once, so a run pinned to
an archived release never touches `rest.ensembl.org`; only if it fails is
`rest.ensembl.org` asked which release is current. Raises `ensembl.ReleaseNotServed` (a
`LookupError`) for a release later than the current one, for an archive reporting another
release, for an archive whose alias redirects to Ensembl's page on archives
(`www.ensembl.org/help/articles/archives`, or the same path on `ensembl.org`:
`ensembl.RETIRED_HOSTS`, `ensembl.RETIRED_PATH`) -- which is how Ensembl retires one;
releases 90 to 104 but 94 did on 2026-09-24, and 100, 103 and 104 on 2026-10-02 -- and
for an archive still unreachable once the retries are spent, an outage or a retirement
that looks like one. A redirect that ends anywhere else off the REST service, such as a
maintenance page, is retried like an outage, not taken for a retirement. A
release after 116 is refused as not on the REST API at all: Ensembl 116 (June 2026) is the
last release of the legacy platform, and its REST API "remains available for e116 for long
term use" with "no plans to port over" to the new one (`ensembl.REST_LAST_RELEASE`). On the
CLI it is `--ensembl-release N` on `annotate` and `identifiability`.

```python
server = ensembl.resolve_server(110)   # 'https://jul2023.rest.ensembl.org' -- GENCODE 44
cfg = annotate.build_config("LEPR", release=110)
res = identifiability.analyze(cfg, ensembl_release=110)
res["annotation"]                      # {'ensembl_release': 110, 'fetched_release': 110,
                                       #  'file_release': None, 'source': None}
```

**Retries.** Every Ensembl request (`annotate`'s lookup included) goes through
`ensembl.request`, or `ensembl.request_json`, which decodes the answer inside the retry
loop. It is retried on HTTP 429/500/502/503/504, connection errors, read timeouts, the
failures below HTTP that `urllib` does not wrap in `URLError` (a connection closed before
the status line, a garbled status line, a body cut short, a reset mid-body) and, for a
JSON request, a 200 whose body is not JSON (`ensembl.BadResponse`); not on any other HTTP
status. `ensembl.TRANSIENT` is the tuple of exceptions that means "the network failed",
for a caller to catch once the retries are spent. `retry` is `retries=` (default
`ensembl.DEFAULT_RETRIES = 5`, counted after the first attempt), `retry_wait=` (default
`ensembl.DEFAULT_RETRY_WAIT = 1.0` s before the first retry, doubling each time; a 429
waits for its `Retry-After` instead) and `timeout=` (default `ensembl.DEFAULT_TIMEOUT =
30` s per attempt). `annotate.build_config`/`run` and `identifiability.analyze` take the
same `retries`/`retry_wait`; on the CLI they are `--retries` and `--retry-wait`. `server=`
selects the host, as returned by `ensembl.resolve_server`.

**`identifiability.kmers(seq, k, canonical=True) -> set`**
Return the set of length-`k` substrings (k-mers) of `seq` (upper-cased), each folded to the
smaller of itself and its reverse complement unless `canonical=False` (the strand-aware
behaviour up to v2.1.1).

**`identifiability.analyze(config, k=31, sequences=None, *, canonical=True, window=None, background_sequences=None, background_fasta=None, decoys=None, background_gene_transcripts="auto", species=None, read_length=100, frag_mean=200.0, frag_sd=60.0, paired=True, depth=30_000_000, mean_efflen=1500.0, tpm=10.0, n_donors=1, conditioning_tau=10.0, min_informative_reads=50.0, min_log2fc=None, ensembl_release=None, inputs_out=None, keep_duplicates=False, max_window_records=20, gtf=None, transcripts_fasta=None, retries=5, retry_wait=1.0) -> dict`**
The three layers the `identifiability` command reports -- sequence uniqueness, the read
model, and estimability on the class-collapsed compatibility system -- for each group and
for the contrast named in `primary_comparison`. `sequences` and `background_sequences`
are optional `{transcript_id: cdna}` dicts; whatever is missing is fetched from Ensembl
release `ensembl_release` (None: the current one), and nothing is fetched when everything
was supplied. Pass a dict as `inputs_out` to get back the sequence the run used
(`sequences`, `background_sequences`, `fetched_release`, `file_release`, `annotation_source`) and the system it was built at
(`k`, `window`, `canonical`, none of which is in the config), which `io.save_inputs`
writes to a file. Returns a dict with `k`, `window`, `canonical`; `annotation`
(`ensembl_release`, `fetched_release`, and `file_release` and `source`, the local files', None without them); `background` (with `gene_id`, the gene the background was fetched as: by the config's `gene_id` when it has one, else by symbol, and a symbol whose gene holds none of the configured transcripts raises `ValueError`; `gene_transcripts` and `n_background_transcripts`, the gene background's columns of the compatibility system; `fasta` and `fasta_sha256`; `fasta_competitors`, the `background_fasta` records that are columns too, each with the number of distinct configured windows it holds, and `n_fasta_competitors`; `max_window_records`; `fasta_left_out`, the records that hold a window of the system and were left out; `windows_dropped`, the distinct windows dropped (`total`, `configured` of `configured_of`, and `columns`, per column id); `sequence_from_fasta`, the gene-background transcripts the FASTA holds with other sequence, whose column is the FASTA's; `decoys`, `decoys_sha256`, `decoys_listed`, `decoys_skipped` and `decoys_absent` (the listed names the FASTA does not hold), all None without `decoys`; `fasta_long_records`, the records read that are longer than `LONG_RECORD` (1 Mb), as id -> length; `same_name_copies`, from `index_scope.fasta_copies`, when `background_fasta` is given, each placed by `index_scope.placed_by_gtf` when `gtf` is; `fasta_other_versions`, the configured transcripts `background_fasta` holds at another version than the sequence used, as `{id: {"used", "fasta"}}` (None without a FASTA); `keep_duplicates`; `identical_to_configured` and `identical_to_background`, the background sequences -- the gene's, `background_sequences`' or FASTA records -- left out because their sequence is a configured transcript's or an earlier background sequence's, as id -> that transcript, since Salmon's default index keeps one of identical sequences, and `identical_source`, `"gene"`, `"sequences"` or `"fasta"` for each -- all three None when `keep_duplicates=True`, which counts every copy); `design`; `groups` (per group:
`verdict`, `reasons`, `n_unique_kmers`, `unique_length`, `n_blocks`,
`expected_informative_reads`, `estimable`, `conditioning_factor`, `gls_relative_se`,
`min_resolvable_log2fc`, `coherence`, ...); `contrast` (the same estimability fields for
the class contrast, plus the effective-length and distinguishing-window summaries);
`counting_noise`; `gene_total`; `effect_resolvable` (None when no `min_log2fc` was
given); `n_compatibility_classes`; `verdict` and `reasons`. The docstring defines each.

A `background_fasta` record with the id of a gene-background transcript is that transcript; every other record that is no copy is an outside record. An outside record that holds a configured window found in at most `max_window_records` outside records (each counted once; default `identifiability.DEFAULT_MAX_WINDOW_RECORDS = 20`) is a column, unless it is longer than `LONG_RECORD`; every window of a column that an outside record left out holds is dropped from every layer -- the unique windows, the read model and the system -- in two passes over the FASTA (the columns' windows, then the added records'). `compatibility_matrix(tracks, ids, drop=)` leaves those windows out and keeps each column's window count as its denominator, so the system is the one with every outside record a column, less the rows that touch a left-out record: what is estimable here is estimable there, with a GLS standard error no smaller. `gene_total` leaves out a column whose every window was dropped and lists it as `transcripts_all_windows_dropped`; `transcripts_without_windows` lists only the transcripts shorter than the window, and it alone sets the CLI's exit 2 (`cli._identifiability_exit`): with windows dropped a column sums to less than one, and the gene total can leave the row space with no transcript windowless. The contrast's `gls_relative_se` is the delta-method SE of log(A/B) only when its gradient `c_a/a - c_b/b` is estimable, and infinite otherwise.

`gtf` and `transcripts_fasta` take the gene's transcripts and their sequence from local files instead of REST, `sequences` and `background_sequences` (`ValueError` with either; on the CLI `--gtf` and `--transcripts-fasta`): GENCODE's comprehensive GTF and the `gencode.vN.transcripts.fa.gz` of the same release, or Ensembl's GTF and cDNA FASTA. The gene is the config's `gene_id`; without one, the one GTF gene that holds every configured transcript; failing that, the symbol's by `annotate.choose_among`. It must hold every configured transcript, and every transcript of it must be in the FASTA at the GTF's version (`annotation_files.AnnotationFileError` otherwise, naming what is missing or which versions differ). The gene background is every other transcript of the gene and goes the way a fetched one does: the identical-copy rule, the id merge with `background_fasta`, `max_window_records` and `decoys`. `transcripts_fasta` alone serves a run with `background_gene_transcripts=False`; it may be the `background_fasta` too. `ensembl_release`, given with `gtf`, must be the release the GTF's header names. Nothing is fetched.

`decoys` is the path of Salmon's `decoys.txt` (one record name per line, read by `io.read_decoys`), for a `background_fasta` that is a decoy-aware index's gentrome: those records are skipped unread. It needs a `background_fasta`; without one it is a `ValueError`.

**`identifiability.scan_fasta_competitors(path, query_kmers, k, canonical=True, exclude_ids=(), identical=None, identical_out=None, *, keep_ids=(), keep_sequences=(), seed=16, decoys=(), stats=None) -> dict`**
The records of a background FASTA that share a window with `query_kmers`, as
`{record id: (sequence, windows)}` in file order; a record that shares none is left out
unless its id is in `keep_ids` or its sequence in `keep_sequences`. `exclude_ids` are
skipped, as are records whose sequence is a key of `identical` (each goes to
`identical_out` as record id -> `identical[sequence]`) and records named in `decoys`.
`stats` receives `read` (False when the query is empty and the file was not read),
`decoys_skipped`, `decoys_found` and `long_records`. Every record reported, among the
competitors or in `identical_out`, has a name of its own: its id; `record<N>` when its
header gives none, N its place among the file's headers; `<id>#2` (`#3`, ...) for a later
record with an id already reported and other sequence (one with the same sequence is the
same record again, and is left out). A record with no sequence is no record.
`identifiability.scan_background_fasta(path, query_kmers, k, ...)` returns the union of
the windows, as through 2.4.1. Each query window is looked up in both orientations, and a
record's window is read only where a `seed`-long substring at a stride of
`k - seed + 1` matches: the answer is the one reading every window gives, and `seed`
changes only the speed.

**`identifiability.estimability(A, c, rcond=1e-10, *, svd=None) -> dict`**
Whether `c'theta` is estimable from `E[y] = A theta` (`estimable`, `residual`, `rank`),
and the structural `conditioning_factor`, from one singular value decomposition with one
tolerance. `svd`, `np.linalg.svd(A.T, full_matrices=False)[:2]`, lets a caller asking about
several functionals of one `A` decompose it once; the answer is the same.

**`identifiability.class_efflen_ratio(lengths_a, lengths_b, frag_mean=200.0) -> (mean_a, mean_b, log2_ratio)`**
Plain mean over each class's transcripts of `effective_length(L, frag_mean) = max(1, L -
frag_mean + 1)`, and the log2 ratio of the two means. Transcript count does not enter; the
docstring records the four alternatives that were scored against simulated bias and why
each lost.

**`identifiability.window_positions(flags) -> list`** and **`position_summary(positions) -> dict`**
Start positions of the flagged windows as fractions of the transcript's own length (0 = 5'
end, 1 = 3' end), and `{n, q1, median, q3}` of them pooled over a class.

The report's `contrast` carries `class_mean_efflen` (`{class: mean}`), `log2_efflen_ratio`
(first class of `primary_comparison` over the second), `efflen_direction_in_band` (whether
`|log2_efflen_ratio| >= DIRECTION_MIN_ABS_LOG2_EFFLEN_RATIO = 0.3`, the band in which a skew
direction was measured) and `distinguishing_window_position` (`{class: position_summary}`
over the windows counted in `n_unique_kmers`).

The full report also carries `annotation`: `{"ensembl_release": <from the config, or
None>, "fetched_release": <the release sequence was fetched from in this run, or None>}`.
`analyze(..., ensembl_release=N)` fetches from release `N` instead of the current one; it
resolves nothing and makes no request when every sequence was supplied.

A group with zero unique k-mers is flagged `distinguishable: False`, which does **not**
mean it cannot be measured: a class with no k-mer of its own is still estimable when a class
it is nested in has unique sequence, so read `verdict` and the estimability fields. Neither
the flag nor `verdict` sets the CLI's exit status: see the exit codes in the README.

---

## `extract`

**`extract.extract(config, quantdir, samplemap, cohort) -> list`**
Read `quantdir/<donor>/quant.sf` Salmon outputs and sum TPM per isoform group.
Returns `[(donor, condition, {group: tpm}), ...]`. Raises `FileNotFoundError` if no
`quant.sf` is found and `ValueError` if a file lacks the `Name`/`TPM` columns.

**`extract.run(config, quantdir, samplemap, cohort, out, allow_mixed_index=False, notes=None) -> int`**
Run `extract` and write a per-donor CSV to `out`; returns the number of donors. Beside it,
`<out>.index.json` (`"format": "isoform-dominance/index/1"`) records each donor's index
from Salmon's `aux_info/meta_info.json` (`index_seq_hash`, `index_name_hash`,
`num_valid_targets`, `keep_duplicates`, `salmon_version`; null for a donor without one),
the distinct `index_seq_hashes`, `mixed`, `missing_meta_info`, `missing_transcripts`
(`n_configured`, `n_missing`, `ids`, `warning`: configured transcripts absent from some
donor's `quant.sf`) and `same_name_copies`. Donors quantified against different indexes
raise `index_scope.MixedIndexError`, and a cohort in which no donor's `quant.sf` has any
configured transcript `index_scope.NoConfiguredTranscripts` (both `io.InputError`s), before
anything is written; the first unless `allow_mixed_index`. Warnings are appended to `notes`
when a list is given.

---

## `annotation_files`

A GTF and a transcript FASTA read offline, with the standard library alone (issue #13).
Ensembl 116 is the last release of the legacy platform and its REST API is kept for 116
only, so an index built from GENCODE 51 or later can be matched only by its own files.

**`annotation_files.scan(path, symbol=None, gene_id=None, transcript_ids=(), block=BLOCK, info=None) -> list`**
Every gene of the GTF named `symbol`, with `gene_id`, or holding one of `transcript_ids`,
as records shaped like REST's `lookup/...?expand=1` -- exactly the fields the package
reads -- so that `annotate.transcripts_of` runs unchanged on either. The file is read in
pieces of `block` bytes (`annotation_files.BLOCK`, 16 MiB), each ending at a line's end, and
only a piece that holds a needle is split into lines. Ids lose their version; a region is
the seqname without `chr`, `chrM` is `MT`; the `_PAR_Y` copies of GENCODE 25-43 are kept
apart (`_par_y`). The protein length is `protein_length`'s; the canonical transcript is the
one tagged exactly `Ensembl_canonical`. A `symbol` no gene has exactly is looked for again
ignoring case (`info["case_insensitive"]`). A GENCODE basic GTF -- every line below a gene
tagged `basic`, in a file of at least `MIN_BASIC_TRANSCRIPTS` (100) transcripts -- raises
`annotation_files.BasicGTF`.

**`annotation_files.protein_length(cds, strand) -> int | None`** — Ensembl's `Translation.length` from a transcript's CDS features `[(start, end, frame)]`, which leave out the stop codon: `(sum of nt + (3 - frame of the 5'-most CDS) % 3) // 3`, the 5'-most CDS on the minus strand being the one that ends furthest along.

**`annotation_files.header(path) -> dict`** — `provider`, `description`, `date`, `gencode_release` (an integer, or `"M37"`-style for mouse) and `ensembl_release`, from GENCODE's `##` header; None for what it does not say (an Ensembl GTF's `#!` header names no release).

**`annotation_files.read_fasta(path, ids) -> dict`** and **`sequences_for(gene, fasta, versions=None) -> dict`** — a transcript FASTA's records by unversioned id, from GENCODE's `|` headers or Ensembl's space-separated ones; `_PAR_Y` records are skipped and an id held twice is an `AnnotationFileError`. `sequences_for` gives every transcript of a `scan` record, and raises `AnnotationFileError` for one the FASTA lacks (with what to pass instead of a subset FASTA) or holds at another version than the GTF.

**`annotation_files.provenance(path) -> dict`** — `{file, bytes, sha256}`. A truncated or corrupt gzip file is one `AnnotationFileError`, not a traceback.

---

## `index_scope`

Same-name copies of a gene on scaffolds, patches and alternate loci, which GENCODE ≥ 48
`transcripts.fa.gz` and Ensembl `cdna.all` contain (README, "Index scope").

**`index_scope.parse_header(head) -> dict | None`**
A GENCODE (`transcript|gene|…|gene name|…`) or Ensembl cDNA (`gene:`, `gene_symbol:`,
`chromosome:`/`scaffold:` region) FASTA header, or a quant.sf name, as
`{transcript, gene_id, gene_name, region}` with unversioned ids (`region` is None for
GENCODE); None for a header in neither format.

**`index_scope.same_name_copies(records, target_ids, gene_names=()) -> list`**
The gene ids among `records` whose gene name is the target transcripts' but whose gene id
is not, as `[{gene_id, gene_name, region, transcripts}]`; empty when no target
transcript is among the records. A record on a reference chromosome
(`index_scope.REFERENCE_REGIONS`: 1–22, X, Y, MT) is left out; only an Ensembl header names
the region, so from a GENCODE header every same-name gene id is returned, with `region`
None, as a possible copy.

**`index_scope.fasta_copies(path, target_ids, gene_names=()) -> list`**
`same_name_copies` over a plain or gzipped FASTA, streamed.

**`index_scope.fasta_copies(..., versions=None)`** also fills `versions`, a dict, with the versioned ids the headers give each target transcript, `{id: {versioned id, ...}}`.

**`index_scope.placed_by_gtf(copies, regions) -> list`** — `copies` with each gene placed by a GTF: `regions` maps the gene ids of the GTF's genes of the name to their region. Each gains `kind`: `"reference_gene"` (on a reference chromosome: the chrY copy of a pseudoautosomal gene, or another gene of the name, which an index of the reference chromosomes keeps), `"off_reference"`, or `"not_in_gtf"` (GENCODE's comprehensive GTF holds the reference chromosomes only, so a copy off them, or a gene of another release); `region` is the GTF's when it holds the gene.

**`index_scope.copy_warning(copies, source, gtf=None) -> str | None`** — the warning the CLI prints;
for copies without a region it says they may be same-name genes on a reference chromosome.
A copy `placed_by_gtf` placed on a reference chromosome is not warned of, and one the GTF
(`gtf`, its name for the message) does not hold is said as such.

**`index_scope.without_identical(copies, identical) -> list`** — `copies` less the records `identical` names (record id -> the transcript it equals: the CLI passes the report's `identical_to_configured` and `identical_to_background`); the CLI warns only about what is left, since Salmon's default index keeps one of identical sequences.

---

## `stats`

**`stats.paired_stat(A, B) -> (n, n_A>B, P, median_fold)`**
Donor-level two-sided Wilcoxon signed-rank test (`scipy.stats.wilcoxon`, `method="auto"`:
exact only without zero differences or ties, otherwise a permutation test at n ≤ 13 and the
normal approximation above; `paired_stat_detail` names which in `wilcoxon_method`) plus
median fold-change. Handles edge cases: empty input → NaNs; all-tied pairs → P is
NaN (test undefined); non-finite ratios are dropped from the fold-change.

**`stats.run(config, condition, cohorts, out) -> dict`**
`cohorts` is `{name: perdonor.csv}`. Writes `<out>.{png,pdf,svg}` and
`<out>_stats.csv`, and returns `{"per_cohort": [...], "combined": (n, n_gt, P, fold),
"detail": [per-cohort paired_stat_detail], "pooled": paired_stat_detail of all donors,
"combination": {"stouffer", "stratified_signed_rank", "pooled"}, "headline_combination":
"stouffer"}`.

---

## `contamination`

**`contamination.run(config, markers, targets, out) -> list`**
Spearman correlation between a contaminant-marker score and the target isoform's
TPM, per cohort (a control for whether a dominance signal is a cell-type artefact).
`markers`/`targets` are `{cohort: path}`. Writes `<out>.{png,pdf,svg}` and
`<out>_scores.csv`; returns rows `(cohort, n, rho, P, median_contam_tissue_ratio)`.

---

## `io`

**`io.load_config(path, need_groups=True) -> dict`** — load a config JSON; `io.InputError` (a `ValueError`) if it is not JSON, not an object, or, unless `need_groups` is false, has no `groups` object.
**`io.load_json(path, what) -> object`** — parse a JSON file; `io.InputError` naming `what` and the file if it is not JSON.
**`io.save_inputs(path, captured, config, release, version, background_fasta=None, decoys=None, max_window_records=None, annotation_source=None) -> dict`** — write the sequence an `identifiability` run used (`captured`, from `analyze(..., inputs_out=...)`) with the release it came from, as `"format": "isoform-dominance/inputs/1"` (`io.INPUTS_FORMAT`). `analysis` records the `k`, `window`, `canonical` and `gene_background` of the run, `max_window_records` when it had a FASTA, and `keep_duplicates` when it had a background, gene or FASTA, none of which the config holds; `sequence_sources` gives each id's `"supplied"`, `"fetched:<release>"` or `"file:<release>"`, and `release` is recorded only when every sequence came from it. A run on local files records them as `annotation_source` (the report's `annotation.source`). A `background_fasta` and a `decoys` file are recorded by path, size and SHA-256, not copied. The gene background is saved whole, copies of configured transcripts included, so that a rerun with `--keep-duplicates` can count them. `gene_id` is the gene the background was fetched as, or the config's `gene_id`. Backs `--save-inputs`.
**`io.load_inputs(path) -> dict`** — read a file `save_inputs` wrote; `io.InputError` if it is not one, and for every field the caller goes on to read: `sequences` and `background_sequences` as objects of id to sequence, an `ensembl_release` that is a number or null, an `analysis` with integer `k` and `window` and boolean `canonical` (and `keep_duplicates`, when present), a `gene_id` that is a string or null, a `background_fasta` and `decoys` that are null or have a path and a sha256, and an `annotation_source` that is null or an object. Which transcripts are a class and which are background is decided by the config the rerun is given, not by the saved grouping. Backs `--inputs`.
**`io.file_sha256(path) -> str`** — hex SHA-256 of a file's bytes.
**`io.read_decoys(path) -> list`** — the record names in Salmon's `decoys.txt`, one per line (its first word, a leading `>` dropped); `io.InputError` when it names none.
**`io.open_text(path)`** — open a text file for reading, gunzipping it when its first two bytes are the gzip magic (`1f 8b`, `io.is_gzip`), whatever its name; `io.open_bytes(path)` the same for bytes.
**`io.shared_transcripts(groups) -> dict`** — `{transcript: [group, ...]}` for each transcript more than one group names; `identifiability` refuses such a config and `extract` warns.
**`io.transcript_to_group(groups) -> dict`** — invert `{group: [ENST...]}` to `{ENST(no version): group}`.
**`io.load_sample_map(path) -> dict`** — read a `donor,condition[,SRR]` CSV to `{donor: condition}`.
**`io.primary_pair(config) -> (gA, gB)`** — the two groups named in `primary_comparison`.
