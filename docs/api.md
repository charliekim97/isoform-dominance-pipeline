# API reference

`isoform-dominance` is primarily a command-line tool, but every subcommand is a
thin wrapper over a small public Python API. Each module is importable from the
`isoform_dominance` package and can be used directly in scripts or notebooks.

```python
from isoform_dominance import (annotate, annotation_files, contamination, ensembl, extract,
                               identifiability, index_scope, io, stats)
```

Every public function and class of those modules is listed below under its module, with
its signature as `inspect.signature` gives it (`tests/test_api_docs.py` holds the page to
the code). A name a module imports from another is listed under the module that defines it.

The config object passed throughout is a plain dict (loaded from JSON via
`io.load_config`) with at least:

```json
{
  "gene": "LEPR",
  "groups": {"short": ["ENST00000371060", "ENST00000616738"], "long": ["ENST00000349533"]},
  "primary_comparison": ["short", "long"]
}
```

A config names one group or more; each lists one or more transcript ids, and
`primary_comparison`, when present, is a list of group names (`io.groups_shape_error`). `annotate` also writes `ensembl_release`,
the Ensembl release the groups were proposed from. It is optional, so configs written
before 2.3.0 still load, but `identifiability` reports a verdict without it as not
reproducible.

---

## `annotate`

Propose isoform groups for a gene from Ensembl by clustering protein-coding
transcripts on their 3' terminal-exon splice acceptor.

**`annotate.AmbiguousGene`** (a subclass of `ValueError`)
A symbol names more than one gene on the reference chromosomes, and none is chosen. The
message lists each gene's id, location and transcript count; `gene_id=` (`--gene-id`)
picks one.

**`annotate.NotOnReference`** (a subclass of `ValueError`)
Raised for a human symbol none of whose genes lies on a reference chromosome (1-22, X, Y,
MT); the message lists each gene with its region, and `gene_id=` takes one of them anyway.
The region and chrX/chrY rules apply to `homo_sapiens` only (`annotate.REFERENCE_SPECIES`);
for any other species two genes of the name are `AmbiguousGene`.

**`annotate.fetch_transcripts(gene, species='homo_sapiens', retries=5, retry_wait=1.0, server=None)`**
Query Ensembl for the gene `lookup/symbol` gives and return `{gene, gene_id, species, strand, transcripts: [{id, protein_aa, terminal_acceptor, is_canonical}]}`. Raises `ValueError` if no protein-coding transcripts with a translation are found. `server` is a base URL from `ensembl.resolve_server`; None is the current release. `retries` and `retry_wait` are as for every request (see Retries under `ensembl`). Network access required.

**`annotate.transcripts_of(g, gene, species='homo_sapiens')`**
What `fetch_transcripts` returns, for an expanded `lookup` record `g`: from REST, or from
a GTF through `annotation_files.scan`, so both sources go through the same code from here.

**`annotate.choose_gene(gene, species, looked_up, **net)`**
The gene a symbol means. `looked_up` is the expanded record `lookup/symbol` gave; `xrefs/symbol` lists the other genes of the name and one `lookup/id` keeps those whose display name is the symbol and that lie on 1–22, X, Y or MT (`xrefs/symbol/SMN1` also lists SMN2 and alternate-locus copies). One such gene: `(it, None)`, with no `lookup/id` when xrefs lists no other gene. A chrX/chrY pair — a pseudoautosomal gene — gives the chrX gene, the copy a GENCODE-built Salmon index keeps (`annotate.GENE_RULE`). Any other two or more raise `annotate.AmbiguousGene` (a `ValueError`) listing each gene's id, location and transcript count. Returns `(record, choice)`: `choice` is `{rule, chosen, candidates, reason}` when there was a choice to make, or when xrefs listed no gene and the others were not looked for.

**`annotate.choose_among(gene, named, species='homo_sapiens')`**
The rule both sources follow, `choose_gene` from REST and `build_config_from_gtf` from a GTF: among `named`, the expanded records of every gene whose display name is the symbol, the human candidates are those on a reference chromosome; one is taken with no choice recorded, of a chrX/chrY pair the chrX gene, and any other two or more raise `AmbiguousGene`. `NotOnReference` when every one of `named` lies off the reference chromosomes. Returns `(record or None, choice or None)`, the record None when no candidate is left.

**`annotate.cluster_by_terminal_exon(info, tolerance=0)`**
Group the transcripts from `fetch_transcripts` by terminal-exon acceptor coordinate; returns clusters, each `{acceptor, rep_aa, n, canonical, ids}`, listed by content: most transcripts first, then the longer representative protein, then the lower acceptor coordinate. So are a config's `_clusters` and `_proposal.tied_with`; through 2.5.0 clusters with as many transcripts were listed in the order the source gave the transcripts, which REST and a GTF give differently. `tolerance` (bp; `annotate --acceptor-tolerance`, issue #4): at 0, the default, a cluster is one exact coordinate. Above 0 the acceptors are taken in ascending order and a cluster takes every acceptor within `tolerance` bp of its first (lowest) one and no further -- nothing chains, so no cluster spans more than `tolerance` -- and its `acceptor` is that first coordinate; a cluster of more than one acceptor also has `members` (`[{acceptor, offset, ids}]`) and `span`. A negative tolerance is a `ValueError`.

**`annotate.alternative_ties(clusters, groups, primary)`**
The clusters the alternative class was chosen over on the tie-breaks rather than on transcript count: those, other than the two proposed, with as many transcripts as the alternative. Empty when the choice was not a tie. A release that adds one transcript to a tied cluster changes the proposal; across six releases of a 109-gene survey, 26 to 39 of 84 to 100 proposals were ties.

**`annotate.propose_groups(info, tolerance=0)`**
Choose the canonical cluster and the largest alternative cluster, returning `(groups, primary, clusters)`: `groups` ({label: [ids]}, labels like `iso_896aa`), the `primary_comparison` list (alternative first), and all clusters. The rule is `annotate.ALTERNATIVE_RULE`: the non-canonical cluster with the most transcripts, a tie going to the longer representative protein, and a tie on both to the lower terminal-acceptor coordinate, so the proposal is a function of the annotation alone.

**`annotate.build_config(gene, species='homo_sapiens', release=None, gene_id=None, acceptor_tolerance=0, **retry)`**
Convenience wrapper returning a complete, reviewable config dict (including `gene_id`, the Ensembl gene the groups were proposed from, `ensembl_release`, `_proposed` notes, `_clusters`, and `_proposal`: `alternative_rule` and `tied_with`, the clusters from `alternative_ties`, and `_gene_choice` from `choose_gene` when there was a choice). `acceptor_tolerance` is `cluster_by_terminal_exon`'s `tolerance`; at 0 the config is the one 2.6.0 wrote, byte for byte. Above 0 the config records `_proposal.acceptor_tolerance` and says it in `_proposed`, and a cluster of more than one acceptor gains `merged_acceptors` (`[{terminal_acceptor, offset, transcripts}]`) and `acceptor_span` in `_clusters`. `gene_id` names the gene outright, by `lookup/id`, and must be a gene of the symbol; it is `--gene-id`. `release` is the Ensembl release to propose the groups from; None is the one `rest.ensembl.org` currently serves, and an earlier one is read from Ensembl's REST archive (see `ensembl.resolve_server`). `ensembl_release` records the release the server reported. `retry` is `retries=` and `retry_wait=`.

**`annotate.build_config_from_gtf(gene, gtf, species='homo_sapiens', gene_id=None, release=None, notes=None, block=16777216, acceptor_tolerance=0)`**
`build_config` from a local GTF -- GENCODE's comprehensive `gencode.vN.annotation.gtf.gz`, or an Ensembl GTF -- with no network (issue #13). The GTF is read by `annotation_files.scan` as REST-shaped records and the gene chosen by `choose_among`, so the config is the one `build_config` writes from the release the file is of, but for `annotation_source`: `{kind: "gtf", file, bytes, sha256, provider, gencode_release, ensembl_release, date, description, n_transcripts}`, `n_transcripts` the gene's transcripts of every biotype. `ensembl_release` is the release the header's `##description` names (`... version 50 (Ensembl 116)`); `release`, when given, must be the same, or `ValueError`, and is recorded as the file's when the header names none. `gene_id` must be a gene of the GTF named `gene`. A symbol no gene has exactly is looked for again ignoring case. A GENCODE basic GTF raises `annotation_files.BasicGTF`. `notes`, a list, receives what a REST run cannot see: a header with no release, a symbol found only ignoring case, a gene none of whose transcripts the GTF tags `Ensembl_canonical` (its canonical class then falls back to the longest protein), and a class of the proposed comparison every one of whose transcripts the GTF tags `mRNA_end_NF`. `block` is the size of the pieces the file is read in (`annotation_files.BLOCK`, 16 MiB). `acceptor_tolerance` is `build_config`'s, through the same function, so the two sources propose alike at every tolerance.

**`annotate.run(gene, out, species='homo_sapiens', release=None, gene_id=None, gtf=None, notes=None, acceptor_tolerance=0, **retry)`**
As `build_config`, or `build_config_from_gtf` when `gtf` is given, but also writes the config JSON to `out` (standard JSON, `io.dump_json`), and returns it. Backs the `annotate` CLI subcommand; `release` is `--ensembl-release`, `gene_id` is `--gene-id`, `gtf` is `--gtf`, `acceptor_tolerance` is `--acceptor-tolerance`, and `retry` is `--retries` and `--retry-wait` (unused with `gtf`).

```python
cfg = annotate.build_config("FLT1")
cfg["groups"]               # {'iso_1338aa': [...], 'iso_733aa': [...]}
cfg["primary_comparison"]   # ['iso_733aa', 'iso_1338aa']
```

---

## `annotation_files`

A GTF and a transcript FASTA read offline, with the standard library alone (issue #13).
Ensembl 116 is the last release of the legacy platform and its REST API is kept for 116
only, so an index built from GENCODE 51 or later can be matched only by its own files.

**`annotation_files.AnnotationFileError`** (a subclass of `io.InputError`)
A GTF or FASTA that cannot give what was asked of it: a transcript missing or at another
version, an id held twice, a truncated or corrupt gzip file.

**`annotation_files.BasicGTF`** (a subclass of `annotation_files.AnnotationFileError`)
A GENCODE `basic` GTF, which leaves out transcripts the index holds.

**`annotation_files.provenance(path)`**
`{file, bytes, sha256}`. A truncated or corrupt gzip file is one `AnnotationFileError`, not a traceback.

**`annotation_files.header(path)`**
`provider`, `description`, `date`, `gencode_release` (an integer, or `"M37"`-style for mouse) and `ensembl_release`, from GENCODE's `##` header; None for what it does not say (an Ensembl GTF's `#!` header names no release).

**`annotation_files.region(seqname)`**
A GENCODE seqname as Ensembl names the region: `chr1` -> `1`, `chrM` -> `MT`.

**`annotation_files.scan(path, symbol=None, gene_id=None, transcript_ids=(), block=16777216, info=None)`**
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

**`annotation_files.incomplete_3prime(transcript)`**
Whether a `scan` transcript is tagged `mRNA_end_NF` (`annotation_files.END_NF`): its 3' end
was not found, so its last annotated exon is not known to be its 3' terminal exon, and is
often an internal exon of a complete transcript. Behind `identifiability`'s
`incomplete_3prime_only` and `annotate --gtf`'s NOTE.

**`annotation_files.incomplete_classes(gene, groups, names)`**
Of the classes `names`, those of `groups` (`{name: [transcript id]}`) every one of whose
transcripts in `gene` (a `scan` record) is `incomplete_3prime`, as `{name: [transcript id]}`.

**`annotation_files.protein_length(cds, strand)`**
Ensembl's `Translation.length` from a transcript's CDS features `[(start, end, frame)]`, which leave out the stop codon: `(sum of nt + (3 - frame of the 5'-most CDS) % 3) // 3`, the 5'-most CDS on the minus strand being the one that ends furthest along.

**`annotation_files.fasta_id(head)`**
The versioned id of a FASTA header (without `>`): the first `|`-field of a GENCODE header,
the first word of an Ensembl one; None when it has none.

**`annotation_files.read_fasta(path, ids)`**
A transcript FASTA's records by unversioned id, `{id: (versioned id, sequence)}` for the records of `ids`, from GENCODE's `|` headers or Ensembl's space-separated ones; `_PAR_Y` records are skipped and an id held twice is an `AnnotationFileError`.

**`annotation_files.sequences_for(gene, fasta, versions=None)`**
`{id: sequence}` for every transcript of a `scan` record, and raises `AnnotationFileError` for one the FASTA lacks (with what to pass instead of a subset FASTA) or holds at another version than the GTF. `versions`, a dict, receives each id's versioned id as the FASTA gives it.

**`annotation_files.version(versioned)`**
The version of a versioned id, as an integer: `ENST00000349533.11` -> 11,
`ENST00000381192.10_PAR_Y` -> 10; None without one.

---

## `identifiability`

Check whether isoform classes, and the contrast between two of them, are measurable by
short reads.

**`identifiability.fetch_cdna(transcript_id, **retry)`**
Fetch one transcript's cDNA sequence from Ensembl, as a string; `ValueError` when Ensembl returns none. For more than one, `ensembl.fetch_cdna_batch`. Network access required.

**`identifiability.fetch_gene(gene, species='homo_sapiens', gene_id=None, **retry)`**
`(gene id, every transcript id)` of a gene, all biotypes, for use as background: by
`gene_id` when given (`lookup/id`), else by the symbol (`lookup/symbol`). Unlike
`annotate`, which keeps protein-coding transcripts only, it keeps the non-coding,
retained-intron and NMD transcripts, which Salmon indexes too. Network access required.

**`identifiability.fetch_gene_transcript_ids(gene, species='homo_sapiens', gene_id=None, **retry)`**
Every transcript id of the gene `fetch_gene` fetches.

**`identifiability.revcomp(seq)`**
The reverse complement of `seq`, upper-cased; `N` stays `N`.

**`identifiability.canonical_kmer(km)`**
The lexicographically smaller of a k-mer and its reverse complement: Salmon's index is built
on canonical k-mers, and an unstranded library gives fragments from either strand.

**`identifiability.kmers(seq, k, canonical=True)`**
Return the set of length-`k` substrings (k-mers) of `seq` (upper-cased), each folded to the
smaller of itself and its reverse complement unless `canonical=False` (the strand-aware
behaviour up to v2.1.1).

**`identifiability.kmer_track(seq, k, canonical=True)`**
The k-mers of `seq` as a list, one per start position, folded as `kmers` folds them; empty
when `seq` is shorter than `k`. What the coverage, block, read-model and system
calculations read, as they need the positions.

**`identifiability.scan_background_fasta(path, query_kmers, k, canonical=True, exclude_ids=(), identical=None, identical_out=None, *, seed=16)`**
The subset of `query_kmers` that also occurs in a background FASTA, plain or gzipped,
streamed one record at a time: the union of the windows `scan_fasta_competitors` reports,
as through 2.4.1. `exclude_ids`, `identical`, `identical_out` and `seed` are as there.

**`identifiability.scan_fasta_competitors(path, query_kmers, k, canonical=True, exclude_ids=(), identical=None, identical_out=None, *, keep_ids=(), keep_sequences=(), seed=16, decoys=(), stats=None)`**
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
Each query window is looked up in both orientations, and a
record's window is read only where a `seed`-long substring at a stride of
`k - seed + 1` matches: the answer is the one reading every window gives, and `seed`
changes only the speed.

**`identifiability.coverage_stats(unique_flags, k)`**
Positional statistics of a transcript's class-unique k-mer starts (`unique_flags[i]` true
when the k-mer at `i` is unique), in bases: `{unique_length, n_blocks, max_block_length,
block_lengths}`. A base is unique when some unique k-mer covers it, so `unique_length` is
the size of the union of the spans, not the sum of their lengths, and a block is a maximal
run of unique bases; `block_lengths` holds the 20 longest, longest first.

**`identifiability.informative_fraction(unique_flags, transcript_length, k, read_length=100, frag_mean=200.0, frag_sd=60.0, paired=True)`**
The probability that a fragment of the transcript is unambiguously assignable to its class:
that a sequenced end -- not the unread middle of the fragment -- covers a whole unique
k-mer. Exact over every start position, averaged over a truncated-normal fragment-length
distribution; 0.0 for a transcript shorter than `k`.

**`identifiability.expected_informative_reads(informative_frac, tpm, transcript_length, depth=30000000, mean_efflen=1500.0, frag_mean=200.0)`**
The expected number of unambiguously assignable fragments,
`depth * tpm * effective_length / (1e6 * mean_efflen) * informative_frac`; `mean_efflen` is
the library's TPM-weighted mean effective length, the one term the gene cannot give.

**`identifiability.effective_length(transcript_length, frag_mean=200.0)`**
`max(1, L - frag_mean + 1)`: the start positions a fragment of mean length has.

**`identifiability.class_efflen_ratio(lengths_a, lengths_b, frag_mean=200.0)`**
Returns `(mean_a, mean_b, log2_ratio)`: the plain mean over each class's transcripts of
`effective_length(L, frag_mean) = max(1, L - frag_mean + 1)`, and the log2 ratio of the
two means; None for all three when either class is empty. Transcript count does not enter;
the docstring records the four alternatives that were scored against simulated bias and
why each lost.

**`identifiability.window_positions(flags)`**
Start positions of the flagged windows as fractions of the transcript's own length (0 = 5'
end, 1 = 3' end).

**`identifiability.position_summary(positions)`**
`{n, q1, median, q3}` of `window_positions` pooled over a class; `n` counts positions, so a
window two transcripts carry counts twice.

The report's `contrast` carries `class_mean_efflen` (`{class: mean}`), `log2_efflen_ratio`
(first class of `primary_comparison` over the second), `efflen_direction_in_band` (whether
`|log2_efflen_ratio| >= DIRECTION_MIN_ABS_LOG2_EFFLEN_RATIO = 0.3`, the band in which a skew
direction was measured) and `distinguishing_window_position` (`{class: position_summary}`
over the windows counted in `n_unique_kmers`).

**`identifiability.unique_read_counting_se(n_informative_a, n_informative_b, n_donors=1)`**
The standard error, from Poisson counting alone, of `log2(n_a / n_b)` for the expected
informative-read counts of each class's best single transcript, as `{log2_ratio_se,
min_resolvable_log2fc, defined}`, `min_resolvable_log2fc` being `1.96 * SE /
sqrt(n_donors)`. Both figures are NaN and `defined` false when either count is 0. No floor
on the contrast's `min_resolvable_log2fc`, which can be smaller.

**`identifiability.counting_noise_floor(n_informative_a, n_informative_b, n_donors=1)`**
Deprecated alias of `unique_read_counting_se`, with a `DeprecationWarning`; removed in 3.0.

**`identifiability.class_coherence(tracks, ids)`**
How alike the transcripts of one class are, as the Jaccard index of their window sets over
every pair: `{median_jaccard, min_jaccard, n_pairs}`, None for both figures with fewer than
two members. The CLI's NOTE on an incoherent class reads `min_jaccard` below 0.05.

**`identifiability.gls_covariance(A, theta)`**
`(A'WA)^+` with `W = diag(1 / A theta)`: the covariance of the GLS estimate of `theta`
under Poisson counts, rows of zero mean dropped; None for an empty `A` or one with no row
of positive mean. `theta` is the expected fragment count of each transcript; `analyze`
passes each transcript's expected fragments at the stated design, the same TPM for all.
Scaling `theta` by `s` for every transcript scales the covariance by `s` and
`sqrt(c'(A'WA)^+ c)` by `sqrt(s)`, so that factor is not compared by size with the
homoskedastic `conditioning_factor`; the docstring gives what does not depend on the scale.

**`identifiability.gls_relative_se(cov, c, theta)`**
The standard error of `c'theta` relative to it, `sqrt(c'cov c) / |c'theta|` -- the SE of
`log(c'theta)`; infinite when `cov` is None, `c` is zero, or the estimate is zero.

**`identifiability.log_ratio_se(cov, c_a, c_b, theta)`**
The delta-method standard error of `log(A/B)` for two class totals from one covariance,
their correlation included; infinite when `cov` is None or either total is zero.

**`identifiability.min_resolvable_log2fc(relative_se, n_donors=1)`**
`1.96 * relative_se / (ln 2 * sqrt(n_donors))`: the |log2 fold change| whose 95% interval
just excludes zero, at about 50% power; infinite when `relative_se` is not finite.

**`identifiability.compatibility_matrix(tracks, transcript_ids, drop=None)`**
The class-collapsed design matrix of the fragment-class system and its classes, `(A,
classes)`: one row per signature (the set of transcripts holding a window), `A[c, t]` the
share of transcript `t`'s window positions in class `c`, and `classes` the sorted ids of
each row's signature. Positions, not distinct windows, are counted, so a column sums to 1
unless the transcript is shorter than the window (all zero) or `drop` removed some of its
windows (the share it kept). A surrogate built from sequence, not the observation model of a
sequencing run: no claim is made relating the two.

**`identifiability.estimability(A, c, rcond=1e-10, *, svd=None)`**
Whether `c'theta` is estimable from `E[y] = A theta` (`estimable`, `residual`, `rank`),
and the structural `conditioning_factor`, from one singular value decomposition with one
tolerance. `svd`, `np.linalg.svd(A.T, full_matrices=False)[:2]`, lets a caller asking about
several functionals of one `A` decompose it once; the answer is the same.

**`identifiability.analyze(config, k=31, sequences=None, *, canonical=True, window=None, background_sequences=None, background_fasta=None, decoys=None, background_gene_transcripts='auto', species=None, read_length=100, frag_mean=200.0, frag_sd=60.0, paired=True, depth=30000000, mean_efflen=1500.0, tpm=10.0, n_donors=1, conditioning_tau=10.0, min_informative_reads=50.0, min_log2fc=None, estimand='fragment', ensembl_release=None, inputs_out=None, keep_duplicates=False, max_window_records=20, gtf=None, transcripts_fasta=None, retries=5, retry_wait=1.0)`**
The three layers the `identifiability` command reports -- sequence uniqueness, the read
model, and estimability on the class-collapsed compatibility system -- for each group and
for the contrast named in `primary_comparison`. A config whose groups or comparison have
the wrong shape is a `ValueError` with the sentence `io.load_config` gives the CLI
(`io.groups_shape_error`), and so is one whose comparison does not name two different
groups of it, or that puts a transcript in two groups. `sequences` and `background_sequences`
are optional `{transcript_id: cdna}` dicts; whatever is missing is fetched from Ensembl
release `ensembl_release` (None: the current one), and nothing is fetched when everything
was supplied. Pass a dict as `inputs_out` to get back the sequence the run used
(`sequences`, `background_sequences`, `fetched_release`, `file_release`, `annotation_source`) and the system it was built at
(`k`, `window`, `canonical`, none of which is in the config), which `io.save_inputs`
writes to a file. Returns a dict with `k`, `window`, `canonical`; `annotation`
(`ensembl_release`, `fetched_release`, and `file_release` and `source`, the local files', None without them); `background` (with `gene_id`, the gene the background was fetched as: by the config's `gene_id` when it has one, else by symbol, and a symbol whose gene holds none of the configured transcripts raises `ValueError`; `gene_transcripts` and `n_background_transcripts`, the gene background's columns of the compatibility system; `fasta` and `fasta_sha256`; `fasta_competitors`, the `background_fasta` records that are columns too, each with the number of distinct configured windows it holds, and `n_fasta_competitors`; `max_window_records`; `fasta_left_out`, the records that hold a window of the system and were left out; `windows_dropped`, the distinct windows dropped (`total`, `configured` of `configured_of`, and `columns`, per column id); `sequence_from_fasta`, the gene-background transcripts the FASTA holds with other sequence, whose column is the FASTA's; `decoys`, `decoys_sha256`, `decoys_listed`, `decoys_skipped` and `decoys_absent` (the listed names the FASTA does not hold), all None without `decoys`; `fasta_long_records`, the records read that are longer than `LONG_RECORD` (1 Mb), as id -> length; `same_name_copies`, from `index_scope.fasta_copies`, when `background_fasta` is given, each placed by `index_scope.placed_by_gtf` when `gtf` is; `fasta_other_versions`, the configured transcripts `background_fasta` holds at another version than the sequence used, as `{id: {"used", "fasta"}}` (None without a FASTA); `keep_duplicates`; `identical_to_configured` and `identical_to_background`, the background sequences -- the gene's, `background_sequences`' or FASTA records -- left out because their sequence is a configured transcript's or an earlier background sequence's, as id -> that transcript, since Salmon's default index keeps one of identical sequences, and `identical_source`, `"gene"`, `"sequences"` or `"fasta"` for each -- all three None when `keep_duplicates=True`, which counts every copy); `design`; `groups` (per group:
`verdict`, `reasons`, `n_unique_kmers`, `unique_length`, `n_blocks`,
`expected_informative_reads`, `estimable`, `conditioning_factor`, `gls_relative_se`,
`finite_se`, `min_resolvable_log2fc`, `coherence`, `incomplete_3prime_only` (since 2.8.0:
with `gtf`, whether the GTF tags every transcript of the class `mRNA_end_NF`, so that the
class may be a fragment of a complete one; None when the transcripts did not come from a
GTF, and it moves no verdict or figure), ...); `contrast` (the same estimability
fields for the class contrast, plus the effective-length and distinguishing-window
summaries); `unique_read_counting_se` (`log2_ratio_se`, `min_resolvable_log2fc`, and
`defined`, false when either class has no informative read and both figures are NaN: the
unique-read counting SE of the log2 class ratio for each class's best transcript, which is
no floor on the contrast's figure; `counting_noise` holds the same, a deprecated key removed
in 3.0, as is the function's old name `counting_noise_floor`); `gene_total`;
`effect_resolvable` (None when no `min_log2fc` was given); `n_compatibility_classes`;
`verdict` and `reasons`; and since 2.9.0, after all of these, `decision_estimand` and
`estimands`. The docstring defines each.

`estimand` (since 2.9.0) is one of `identifiability.ESTIMANDS`, `"fragment"` (the default) or
`"abundance"`, and picks the class totals that decide `verdict`, `reasons`, `effect_resolvable`
and the estimand's fields of `groups` and `contrast` (`estimable`, `residual`,
`conditioning_factor`, `rank`, `label`, `gls_relative_se`, `finite_se`,
`min_resolvable_log2fc`, `beyond_linear`, `verdict`, `reasons`). The fragment estimand is a
class's expected fragment count, the sum of theta_t over its transcripts: 2.8.0's question, and
what a quantifier's `NumReads` sums. The abundance estimand is its length-corrected abundance,
the sum of theta_t / l_t with l_t the transcript's window positions (its effective length in
this system), whose class ratio is the ratio that summed TPM estimates -- what `extract` and
`stats` compare. They are different functionals of one system and can part: a direction the
system cannot see that moves abundance between transcripts of one class that differ in length
can leave the class's abundance where it was and move its fragment count, and then the fragment
estimand is not estimable, while the abundance one is when no direction the system cannot see
changes a compared class's abundance (`identifiability.ESTIMANDS` says more). Its coefficients,
1/l_t, are scaled to the norm of the fragment functional's, which moves no estimability and no
relative standard error, only the conditioning factor, onto the fragment one's footing. Both
are computed and given under `estimands`, each as
`{definition, precision_checked_in_simulation, groups, contrast, effect_resolvable, verdict, reasons}`
with only the fields above per class and for the contrast; `decision_estimand` names the one
that decides. With the default, every field 2.8.0 wrote is the fragment estimand's, as it was.
Where only the abundance estimand is estimable, that answer rests on counting a transcript's
length in window positions, and its figures weight a short transcript less than TPM does
(`identifiability.ESTIMANDS` gives both). Its precision figures have not yet been checked
against a quantifier's TPM in simulation (its `precision_checked_in_simulation` is `False`, the
fragment one's `True`); `validation/README.md` says what that check must show before it becomes
the default in 3.0. Any other value is a `ValueError`.

`min_resolvable_log2fc` is `1.96 * SE / sqrt(n_donors)` on the log2 scale: the |log2 fold
change| whose estimate's 95% interval just excludes zero. A true effect of exactly that
size is detected about half the time -- about 50% power; for 80% power the effect has to
be about `2.8 * SE / sqrt(n_donors)`, 1.43 times the figure. No option computes that.

A figure that is not finite is infinite in the returned dict and null in `--json`, which is
standard JSON (`io.finite_json`). `estimable` false says why a `conditioning_factor` is;
`finite_se` false says why `gls_relative_se` and `min_resolvable_log2fc` are, including for
a contrast that is estimable while its class totals are not, which `estimable` cannot say;
`unique_read_counting_se.defined` says it for the unique-read counting SE.

A `background_fasta` record with the id of a gene-background transcript is that transcript; every other record that is no copy is an outside record. An outside record that holds a configured window found in at most `max_window_records` outside records (each counted once; default `identifiability.DEFAULT_MAX_WINDOW_RECORDS = 20`) is a column, unless it is longer than `LONG_RECORD`; every window of a column that an outside record left out holds is dropped from every layer -- the unique windows, the read model and the system -- in two passes over the FASTA (the columns' windows, then the added records'). `compatibility_matrix(tracks, ids, drop=)` leaves those windows out and keeps each column's window count as its denominator, so the system is the one with every outside record a column, less the rows that touch a left-out record: what is estimable here is estimable there, with a GLS standard error no smaller. `gene_total` leaves out a column whose every window was dropped and lists it as `transcripts_all_windows_dropped`; `transcripts_without_windows` lists only the transcripts shorter than the window, and it alone sets the CLI's exit 2 (`cli._identifiability_exit`): with windows dropped a column sums to less than one, and the gene total can leave the row space with no transcript windowless. The contrast's `gls_relative_se` is the delta-method SE of log(A/B) only when its gradient `c_a/a - c_b/b` is estimable, and infinite otherwise (`finite_se` false; null in `--json`).

`gtf` and `transcripts_fasta` take the gene's transcripts and their sequence from local files instead of REST, `sequences` and `background_sequences` (`ValueError` with either; on the CLI `--gtf` and `--transcripts-fasta`): GENCODE's comprehensive GTF and the `gencode.vN.transcripts.fa.gz` of the same release, or Ensembl's GTF and cDNA FASTA. The gene is the config's `gene_id`; without one, the one GTF gene that holds every configured transcript; failing that, the symbol's by `annotate.choose_among`. It must hold every configured transcript, and every transcript of it must be in the FASTA at the GTF's version (`annotation_files.AnnotationFileError` otherwise, naming what is missing or which versions differ). The gene background is every other transcript of the gene -- and `background.gene_id` the gene, as for a fetched background -- and goes the way a fetched one does: the identical-copy rule, the id merge with `background_fasta`, `max_window_records` and `decoys`. `transcripts_fasta` alone serves a run with `background_gene_transcripts=False`; it may be the `background_fasta` too. `ensembl_release`, given with `gtf`, must be the release the GTF's header names. Nothing is fetched.

`decoys` is the path of Salmon's `decoys.txt` (one record name per line, read by `io.read_decoys`), for a `background_fasta` that is a decoy-aware index's gentrome: those records are skipped unread. It needs a `background_fasta`; without one it is a `ValueError`.

The full report also carries `annotation`: `{"ensembl_release": <from the config, or
None>, "fetched_release": <the release sequence was fetched from in this run, or None>}`.
`analyze(..., ensembl_release=N)` fetches from release `N` instead of the current one; it
resolves nothing and makes no request when every sequence was supplied.

```python
server = ensembl.resolve_server(110)   # 'https://jul2023.rest.ensembl.org' -- GENCODE 44
cfg = annotate.build_config("LEPR", release=110)
res = identifiability.analyze(cfg, ensembl_release=110)
res["annotation"]                      # {'ensembl_release': 110, 'fetched_release': 110,
                                       #  'file_release': None, 'source': None}
```

A group with zero unique k-mers is flagged `distinguishable: False`, which does **not**
mean it cannot be measured: a class with no k-mer of its own is still estimable when a class
it is nested in has unique sequence, so read `verdict` and the estimability fields. Neither
the flag nor `verdict` sets the CLI's exit status: see the exit codes in the README.

**`identifiability.run(config, k=31, sequences=None, **kw)`**
A thin wrapper of `analyze`, kept for backward compatibility.

---

## `ensembl`

Every request to the Ensembl REST API, retried, and the release it serves.

**`ensembl.BadResponse`** (a subclass of `urllib.error.URLError`)
A 200 whose body is not what the request asked for, such as HTML where JSON was due: a
`URLError`, so it is retried and reported as a network failure, not as a config problem.

**`ensembl.request(path, body=None, *, server=None, accept='application/json', timeout=30, retries=5, retry_wait=1.0)`**
GET `path`, or POST `body` to it as JSON, and return the response bytes. `server` is the
base URL, `ensembl.SERVER` when None; for an earlier release pass what `resolve_server`
returns, never the `e<N>` alias itself. The last error is raised once `retries` retries
have been spent.

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

**`ensembl.request_json(path, body=None, **retry)`**
As `request`, decoding the JSON answer inside the retry loop: a 200 that is not JSON raises
`BadResponse` and is retried like any other transient failure.

**`ensembl.get_json(path, **retry)`**
GET `path` and decode the JSON response; `retry` goes to `request`.

**`ensembl.release_number(info)`**
The release in an `/info/data` response, which lists the releases a server holds;
`ValueError` unless it lists exactly one, because the number is recorded as the release a
verdict was computed on.

**`ensembl.fetch_release(**retry)`**
The release the server is serving, from `GET /info/data`; raises `ValueError` unless it
lists exactly one. `annotate.build_config` records it, and `analyze` reports it as
`annotation.fetched_release` when, and only when, the run fetched sequence.

**`ensembl.ReleaseNotServed`** (a subclass of `LookupError`)
No REST server answers for the Ensembl release that was asked for; see `resolve_server`.

**`ensembl.resolve_server(release=None, *, timeout=30, retries=5, retry_wait=1.0)`**
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
(`https://www.ensembl.org/help/articles/archives`, or the same path on the bare host:
`ensembl.RETIRED_HOSTS`, `ensembl.RETIRED_PATH`) -- which is how Ensembl retires one;
releases 90 to 104 but 94 did on 2026-09-24, and 100, 103 and 104 on 2026-10-02 -- and
for an archive still unreachable once the retries are spent, an outage or a retirement
that looks like one. A redirect that ends anywhere else off the REST service, such as a
maintenance page, is retried like an outage, not taken for a retirement. A
release after 116 is refused as not on the REST API at all: Ensembl 116 (June 2026) is the
last release of the legacy platform, and its REST API "remains available for e116 for long
term use" with "no plans to port over" to the new one (`ensembl.REST_LAST_RELEASE`). On the
CLI it is `--ensembl-release N` on `annotate` and `identifiability`. `timeout`, `retries`
and `retry_wait` are keyword-only, as for `request`.

**`ensembl.fetch_cdna_batch(ids, versions=None, **retry)`**
`{id: cdna}` via `POST /sequence/id`, 50 ids per request, matched back to ids by the
`query` field Ensembl echoes. `analyze` fetches all of its sequence this way. `versions`, a
dict, receives each sequence's versioned id when the answer gives one. An id Ensembl
returns nothing for is absent from the result.

---

## `extract`

**`extract.NoQuantFiles`** (a subclass of `io.InputError` and `FileNotFoundError`)
No `<quantdir>/<donor>/quant.sf`: an `io.InputError` for the CLI's one line, and a
`FileNotFoundError` for callers that caught that before 2.4.1.

**`extract.quant_paths(quantdir)`**
`{donor: quantdir/<donor>/quant.sf}`; `NoQuantFiles` if there is none.

**`extract.extract(config, quantdir, samplemap, cohort, found_out=None)`**
Read `quantdir/<donor>/quant.sf` Salmon outputs and sum TPM per isoform group.
Returns `[(donor, condition, {group: tpm}), ...]`, in the order of the quant.sf paths. Raises
`NoQuantFiles` (a `FileNotFoundError`) if no `quant.sf` is found and `ValueError` if a file
lacks the `Name`/`TPM` columns. `found_out`, a dict, receives each donor's set of configured
transcript ids its `quant.sf` has.

**`extract.missing_transcripts(config, found, cohort)`**
What the sidecar records about configured transcripts absent from some donor's `quant.sf`
(`found` is `extract`'s `found_out`): `{n_configured, n_missing, ids, warning}`, `warning`
one line for the cohort or None. Raises `index_scope.NoConfiguredTranscripts` when no
donor has any of them, as every class total would be zero.

**`extract.write_perdonor(config, rows, cohort, out)`**
Writes `extract`'s rows to the per-donor CSV `out` and returns `out`: `cohort, donor,
condition`, a `<group>_TPM` column per group, and `<A>_fraction`, `A / (A + B)` for the two
groups of `primary_comparison` (without one, the first two groups; `NA` when both are zero),
left out unless the comparison names exactly two different groups of the config;
`extract.run` warns when it names a group the config lacks, other than two groups, or one
group twice.

**`extract.quant_copies(config, prov, quants)`**
Same-name copies among the targets of each index the cohort was quantified against. Only
an index built without Salmon's `--gencode` keeps GENCODE's whole header as the target
name, and only then does `quant.sf` say which gene each target belongs to; plain names give
nothing to check, and nothing is reported.

**`extract.run(config, quantdir, samplemap, cohort, out, allow_mixed_index=False, notes=None)`**
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

## `stats`

**`stats.load_perdonor(path, condition, gA, gB)`**
Reads a per-donor CSV as `extract` writes it: `(donors, A, B)`, the `<gA>_TPM` and
`<gB>_TPM` columns as arrays, of the rows whose `condition` is `condition` (every row when
it is None, `""` or `"all"`); `io.InputError` naming a column the table lacks.

**`stats.signed_rank_resolution_floor(n, alternative='two-sided')`**
The smallest p-value the signed-rank sign-permutation null can give with `n` non-zero
pairs: `2^(1-n)` two-sided, `2^-n` otherwise; NaN for `n` of 0 or less. At n = 5 it is
0.0625, so no five donors reach 0.05. Ties among the absolute differences do not raise it.

**`stats.paired_stat(A, B)`**
Returns `(n, n_A>B, P, median_fold)`: donor-level two-sided Wilcoxon signed-rank test (`scipy.stats.wilcoxon`, `method="auto"`:
exact only without zero differences or ties, otherwise a permutation test at n ≤ 13 and the
normal approximation above; `paired_stat_detail` names which in `wilcoxon_method`) plus
median fold-change. Handles edge cases: empty input → NaNs; all-tied pairs → P is
NaN (test undefined); non-finite ratios are dropped from the fold-change.

**`stats.signed_rank_direction(d, zero_method='wilcox')`**
+1.0, -1.0 or 0.0: the sign of `W+ - W-` with the ranks the test uses (non-zero differences
under `wilcox`; all under `pratt` and `zsplit`), that is the side of its null centre the
statistic that produced the p-value fell on. `d = [1, 1, -5, -6]` is -1.0 (median -2, `W+`
3 of 10), where the count of positive pairs is a tie.

**`stats.paired_stat_detail(A, B, n_boot=10000, seed=0, zero_method='wilcox', ci=0.95)`**
The same test with what a reader needs beside it, as a dict: `n`, `n_greater`, `p`, `wilcoxon_method`
(the computation that produced `p`: `exact`, `permutation`, `asymptotic`, or `unresolved`),
`median_fold`; `n_ties` (zero differences), `n_effective` (the pairs the test ranks: the
non-zero ones under `wilcox`, all under `pratt` and `zsplit`), `n_nonzero`; `direction`
(`signed_rank_direction`); `p_floor`, `2^(1 - n_nonzero)` under every zero method (a zero
has no sign to flip, so `pratt` and `zsplit`, which rank zeros, reach no lower: through
2.6.0 their floor counted every pair), `underpowered` (floor > 0.05) and
`approximation_below_floor` (`p` below the floor, which only a normal approximation can
give); `fold_ci` and `fold_ci_method`: the donor-bootstrap 95% percentile interval of the
median ratio (`"bootstrap"`), or, from five finite ratios or fewer, their minimum and
maximum (`"range (n<=5)"`) -- which is what the bootstrap interval of such a median is,
since the bootstrap median equals the minimum (and the maximum) with probability 0.25,
7/27, 13/256 and 181/3125 at n = 2 to 5, each above 0.025, and 0.0087 at n = 6
(`stats.RANGE_MAX_N`). With `n_boot=0` there is no interval and
`fold_ci_method` is None.

**`stats.stouffer(per_cohort)`**
`per_cohort` is `(name, n, direction, p)` or `(name, n, direction, p, p_floor)` per cohort.
Each two-sided p becomes `direction * norm.isf(p / 2)` -- `p / 2` held at the smallest
normal double, so a p of 0 gives z = 37.5, not inf -- weighted by `sqrt(n)`; `run` passes
`n_effective` and `direction`. Returns `{z, p, k, p_floor, underpowered,
approximation_below_floor}`: `p_floor` is every cohort at its floor in one direction,
`2 sf(sum(w_i z_i) / sqrt(sum(w_i^2)))` with `z_i` each floor's z, and NaN unless every
cohort's floor is given. Through 2.6.0 `run` signed a cohort by its count of positive pairs
(ties to +) and weighted it by `sqrt(n)` of every donor.

**`stats.stratified_signed_rank(strata)`**
Within each stratum of paired differences, `W+` over the non-zero pairs with average ranks,
centred by `n(n+1)/4` and scaled by the tie-corrected null variance
`n(n+1)(2n+1)/24 - sum(t^3 - t)/48` (through 2.6.0 without the tie term); strata weighted
`1/(n+1)` and combined by the normal approximation. Returns `{z, p, k, p_floor,
underpowered, approximation_below_floor, note}`: `p_floor` is `2^(1 - N)`, `N` the
non-zero pairs of all strata, and `note` is `stats.BELOW_FLOOR_NOTE` ("approximation below
the attainable exact minimum") when `p` lies below it -- one stratum `d = 1..5` gives
0.0431 against 0.0625.

**`stats.run(config, condition, cohorts, out, n_boot=10000, seed=0, zero_method='wilcox')`**
`cohorts` is `{name: perdonor.csv}`. Writes `<out>.{png,pdf,svg}` and
`<out>_stats.csv`, and returns `{"per_cohort": [...], "combined": (n, n_gt, P, fold),
"detail": [per-cohort paired_stat_detail], "pooled": paired_stat_detail of all donors,
"combination": {"stouffer", "stratified_signed_rank", "pooled"}, "headline_combination":
"stouffer"}`; each combination carries `p_floor`, `underpowered` and
`approximation_below_floor`, and `pooled`'s `z` is NaN (null in `--json`), as it is no
z-combination. The table has one row per cohort, `COMBINED_POOLED`,
`COMBINED_STOUFFER` and `COMBINED_STRATIFIED_SIGNED_RANK`, with the columns `cohort, n,
A>B, median_fold, fold_CI_low, fold_CI_high, paired_wilcoxon_P, resolution_floor_P,
underpowered, wilcoxon_method, fold_CI_method, note`; every row has its floor, and the last
two columns are new in 2.7.0, after the old ones. The two groups are `io.primary_pair`'s.

---

## `contamination`

**`contamination.load_markers(path, tissue, contaminant)`**
Reads a marker-TPM CSV: `{donor: (tissue score, contaminant score, contaminant / tissue)}`,
each score the mean `log2(TPM + 1)` over the panel's genes the table has (the ratio NaN
when the tissue score is 0); `io.InputError` when it has no gene of either panel.

**`contamination.load_target(path, target_group)`**
`{donor: TPM}` from the `<target_group>_TPM` column of a per-donor CSV; `io.InputError`
when the table has no such column.

**`contamination.qc_section(config)`**
`(target_group, tissue markers, contaminant markers)` from the config's `contamination_qc`;
`io.InputError` saying what is missing.

**`contamination.run(config, markers, targets, out)`**
Spearman correlation between a contaminant-marker score and the target isoform's
TPM, per cohort (a control for whether a dominance signal is a cell-type artefact).
`markers`/`targets` are `{cohort: path}`. Writes `<out>.{png,pdf,svg}` and
`<out>_scores.csv`; returns rows `(cohort, n, rho, P, median_contam_tissue_ratio)`.

---

## `io`

Every JSON the package writes -- each `--json`, the `--save-inputs` file and the config
`annotate` writes -- is standard JSON (RFC 8259): `io.finite_json` returns its argument with
every number that is not finite as None, and `io.dump_json` writes with `allow_nan=False`,
so a NaN or infinity that reaches a writer is a `ValueError`, not a `NaN` or `Infinity`
token a standard parser refuses. Reading is lenient: `io.load_json` reads those tokens as
Python does.

**`io.InputError`** (a subclass of `ValueError`)
A file named on the command line that cannot be used as it is; the CLI says it in one
`error:` line and exits 1.

**`io.load_json(path, what)`**
Parse a JSON file; `io.InputError` naming `what` and the file if it is not JSON, or not UTF-8.

**`io.finite_json(obj)`**
`obj` with every number that is not finite -- NaN, inf, -inf -- as None, through dicts,
lists and tuples (a tuple becomes a list).

**`io.dump_json(obj, path, **kw)`**
Write `obj` to `path` with `allow_nan=False`; `kw` goes to `json.dump`.

**`io.is_gzip(path)`**
Whether the file's first two bytes are the gzip magic (`1f 8b`), whatever its name.

**`io.open_text(path)`**
Open a text file for reading, gunzipping it when `io.is_gzip`, whatever its name.

**`io.open_bytes(path)`**
The same as `open_text`, for bytes.

**`io.load_config(path, need_groups=True)`**
Load a config JSON; `io.InputError` (a `ValueError`) if it is not JSON, not an object, or, unless `need_groups` is false (`qc`), has no `groups` object or fails `io.groups_shape_error`: then the message is `config <path>: ` and that sentence, which names `groups`, the group or `primary_comparison` at fault. `identifiability`, `extract` and `stats` read the config here.

**`io.groups_shape_error(config)`**
What is wrong with the shape of a config's groups and comparison, as a sentence that starts
with the place (`groups`, `group "<name>"` or `primary_comparison`), or None: `groups` must
name one group or more, each group a non-empty list of transcript ids, each a non-empty
string with no white space, and `primary_comparison`, when present, a list of group names,
each a string. From Python a list may be any collection that is not a string or a mapping
(a tuple, a set, an array), but a comparison has an order and so is no set. Since 2.9.0;
through 2.8.0 a group given as a string was read one character at a time as transcript ids,
an empty group gave a verdict, an id with a space was read as it stood (unless a version
followed it), a group given as an object was read as its keys, and a comparison given as a
string was read as its letters. That the comparison names exactly two groups, both the
config's and not one twice, is checked by `io.primary_pair` (for `stats`) and
`identifiability.analyze`; `extract` writes no fraction column for another comparison, and
says so.

**`io.file_sha256(path, chunk=1048576)`**
Hex SHA-256 of a file's bytes, read `chunk` bytes (1 MiB) at a time.

**`io.read_decoys(path)`**
The record names in Salmon's `decoys.txt`, one per line (its first word, a leading `>` dropped); `io.InputError` when it names none.

**`io.transcript_to_group(groups)`**
Invert `{group: [ENST...]}` to `{ENST(no version): group}`; a transcript in two groups goes to the one listed last.

**`io.shared_transcripts(groups)`**
`{transcript: [group, ...]}` for each transcript more than one group names; `identifiability` refuses such a config and `extract` warns.

**`io.load_sample_map(path)`**
Read a `donor,condition[,SRR]` CSV to `{donor: condition}`; `io.InputError` without a `donor` column.

**`io.primary_pair(config)`**
The two groups named in `primary_comparison`, `(gA, gB)`; without one, the config's first
two groups. A `ValueError` when it does not name exactly two groups (through 2.8.0 a third
and later were ignored), names a group the config lacks (through 2.8.0 left to the
per-donor table, which named the missing column) or names one group twice, as
`identifiability` refuses such a comparison, and for a set, which has no order. `stats`
reads the pair here, and the CLI refuses such a config in one `error: config <path>: ...`
line before reading a table.

**`io.save_inputs(path, captured, config, release, version, background_fasta=None, decoys=None, max_window_records=None, annotation_source=None)`**
Write the sequence an `identifiability` run used (`captured`, from `analyze(..., inputs_out=...)`) with the release it came from, as `"format": "isoform-dominance/inputs/1"` (`io.INPUTS_FORMAT`), and return the document. `analysis` records the `k`, `window`, `canonical` and `gene_background` of the run, `max_window_records` when it had a FASTA, and `keep_duplicates` when it had a background, gene or FASTA, none of which the config holds; `sequence_sources` gives each id's `"supplied"`, `"fetched:<release>"` or `"file:<release>"`, and `release` is recorded only when every sequence came from it. A run on local files records them as `annotation_source` (the report's `annotation.source`). `versions` gives each sequence's versioned id where the run knew it, and a rerun passes it on, so `background.fasta_other_versions` survives `--inputs`. A `background_fasta` and a `decoys` file are recorded by path, size and SHA-256, not copied. The gene background is saved whole, copies of configured transcripts included, so that a rerun with `--keep-duplicates` can count them. `gene_id` is the gene the background was fetched as, or the config's `gene_id`. Backs `--save-inputs`.

**`io.load_inputs(path)`**
Read a file `save_inputs` wrote; `io.InputError` if it is not one, and for every field the caller goes on to read: `sequences` and `background_sequences` as objects of id to sequence, an `ensembl_release` that is a number or null, an `analysis` with integer `k` and `window` and boolean `canonical` (and `keep_duplicates`, when present), a `gene_id` that is a string or null, a `background_fasta` and `decoys` that are null or have a path and a sha256, an `annotation_source` that is null or names the files a run read (`kind`, `gtf`, `transcripts_fasta`), and `versions` that is null or an object of id to versioned id. Which transcripts are a class and which are background is decided by the config the rerun is given, not by the saved grouping. Backs `--inputs`.

---

## `index_scope`

Same-name copies of a gene on scaffolds, patches and alternate loci, which GENCODE ≥ 48
`transcripts.fa.gz` and Ensembl `cdna.all` contain (README, "Index scope"), and which
Salmon index quantified each donor.

**`index_scope.MixedIndexError`** (a subclass of `io.InputError`)
Donors in one cohort were quantified against different indexes.

**`index_scope.NoConfiguredTranscripts`** (a subclass of `io.InputError`)
No donor's `quant.sf` has any transcript the config names.

**`index_scope.parse_header(head)`**
A GENCODE (`transcript|gene|…|gene name|…`) or Ensembl cDNA (`gene:`, `gene_symbol:`,
`chromosome:`/`scaffold:` region) FASTA header, or a quant.sf name, as
`{transcript, gene_id, gene_name, region}` with unversioned ids (`region` is None for
GENCODE); None for a header in neither format.

**`index_scope.same_name_copies(records, target_ids, gene_names=())`**
The gene ids among `records` whose gene name is the target transcripts' but whose gene id
is not, as `[{gene_id, gene_name, region, transcripts}]`; empty when no target
transcript is among the records. A record on a reference chromosome
(`index_scope.REFERENCE_REGIONS`: 1–22, X, Y, MT) is left out; only an Ensembl header names
the region, so from a GENCODE header every same-name gene id is returned, with `region`
None, as a possible copy.

**`index_scope.copies_in(read_heads, target_ids, gene_names=(), versions=None)`**
`same_name_copies` over the headers `read_heads()` yields, holding only the records that
can matter rather than every header of a whole-transcriptome FASTA. `read_heads` is called
once more when the name the headers give the target transcripts is not among `gene_names`:
a copy can come before its reference gene. `versions` is as for `fasta_copies`.

**`index_scope.fasta_copies(path, target_ids, gene_names=(), versions=None)`**
`same_name_copies` over a plain or gzipped FASTA, streamed (`copies_in`). `versions`, a dict, is filled with the versioned ids the headers give each target transcript, `{id: {versioned id, ...}}`.

**`index_scope.placed_by_gtf(copies, regions)`**
`copies` with each gene placed by a GTF: `regions` maps the gene ids of the GTF's genes of the name to their region. Each gains `kind`: `"reference_gene"` (on a reference chromosome: the chrY copy of a pseudoautosomal gene, or another gene of the name, which an index of the reference chromosomes keeps), `"off_reference"`, or `"not_in_gtf"` (GENCODE's comprehensive GTF holds the reference chromosomes only, so a copy off them, or a gene of another release); `region` is the GTF's when it holds the gene.

**`index_scope.without_identical(copies, identical)`**
`copies` less the records `identical` names (record id -> the transcript it equals: the CLI passes the report's `identical_to_configured` and `identical_to_background`); the CLI warns only about what is left, since Salmon's default index keeps one of identical sequences.

**`index_scope.copy_warning(copies, source, gtf=None)`**
The warning the CLI prints, or None;
for copies without a region it says they may be same-name genes on a reference chromosome.
A copy `placed_by_gtf` placed on a reference chromosome is not warned of, and one the GTF
(`gtf`, its name for the message) does not hold is said as such.

**`index_scope.UnreadableMetaInfo`** (a subclass of `Exception`)
A `meta_info.json` that is there but cannot be read as Salmon's; the message names it.

**`index_scope.read_meta_info(quant_path)`**
The index fields of the Salmon `aux_info/meta_info.json` beside a `quant.sf`, or None when
there is none (kallisto, or an older Salmon); `UnreadableMetaInfo` when the file is there
but is not a JSON object in UTF-8 whose hashes are strings.

**`index_scope.index_provenance(donor_quants)`**
`{donor: quant.sf path}` -> the sidecar document: each donor's index fields, the distinct
`index_seq_hash` and `index_decoy_seq_hash` values, whether there is more than one index,
the donors with no `meta_info.json`, and those whose `meta_info.json` could not be read
(`{donor: why}`). Two donors are on different indexes when their `index_seq_hash` differs,
or when both record an `index_decoy_seq_hash` and it differs.

**`index_scope.mixed_index_message(prov)`**
Which donors went with which index, from `index_provenance`'s document: what
`MixedIndexError`'s message and the warning `--allow-mixed-index` gives say of it.

---

## `cli`

Use the commands as `isoform-dominance` and its subcommands; the functions of this module
are internal to that command line, not part of the API.
