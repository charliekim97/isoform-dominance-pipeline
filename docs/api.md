# API reference

`isoform-dominance` is primarily a command-line tool, but every subcommand is a
thin wrapper over a small public Python API. Each module is importable from the
`isoform_dominance` package and can be used directly in scripts or notebooks.

```python
from isoform_dominance import annotate, identifiability, extract, stats, contamination, io
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
Query Ensembl and return `{gene, species, strand, transcripts: [{id, protein_aa, terminal_acceptor, is_canonical}]}`. Raises `ValueError` if no protein-coding transcripts with a translation are found. `server` is a base URL from `ensembl.resolve_server`; None is the current release. Network access required.

**`annotate.cluster_by_terminal_exon(info) -> list`**
Group the transcripts from `fetch_transcripts` by terminal-exon acceptor coordinate; returns clusters sorted by size, each `{acceptor, rep_aa, n, canonical, ids}`.

**`annotate.propose_groups(info) -> (groups, primary, clusters)`**
Choose the canonical cluster and the largest alternative cluster, returning `groups` ({label: [ids]}, labels like `iso_896aa`), the `primary_comparison` list (alternative first), and all clusters. The rule is `annotate.ALTERNATIVE_RULE`: the non-canonical cluster with the most transcripts, a tie going to the longer representative protein, and a tie on both to the lower terminal-acceptor coordinate, so the proposal is a function of the annotation alone.

**`annotate.alternative_ties(clusters, groups, primary) -> list`**
The clusters the alternative class was chosen over on the tie-breaks rather than on transcript count: those, other than the two proposed, with as many transcripts as the alternative. Empty when the choice was not a tie. A release that adds one transcript to a tied cluster changes the proposal; across six releases of a 109-gene survey, 26 to 39 of 84 to 100 proposals were ties.

**`annotate.build_config(gene, species="homo_sapiens", release=None) -> dict`**
Convenience wrapper returning a complete, reviewable config dict (including `ensembl_release`, `_proposed` notes, `_clusters`, and `_proposal`: `alternative_rule` and `tied_with`, the clusters from `alternative_ties`). `release` is the Ensembl release to propose the groups from; None is the one `rest.ensembl.org` currently serves, and an earlier one is read from Ensembl's REST archive (see `ensembl.resolve_server`). `ensembl_release` records the release the server reported.

**`annotate.run(gene, out, species="homo_sapiens", release=None) -> dict`**
As `build_config`, but also writes the config JSON to `out`. Backs the `annotate` CLI subcommand; `release` is `--ensembl-release`.

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

**`ensembl.fetch_cdna_batch(ids, **retry) -> dict`**
`{id: cdna}` via `POST /sequence/id`, 50 ids per request, matched back to ids by the
`query` field Ensembl echoes. `analyze` fetches all of its sequence this way.

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
release, for an archive whose alias redirects off the REST service -- which is how Ensembl
retires one; releases 90 to 104 but 94 did on 2026-09-24 -- and for an archive still
unreachable once the retries are spent, an outage or a retirement that looks like one. A
release after 116 is refused as not on the REST API at all: Ensembl 116 (June 2026) is the
last release of the legacy platform, and its REST API "remains available for e116 for long
term use" with "no plans to port over" to the new one (`ensembl.REST_LAST_RELEASE`). On the
CLI it is `--ensembl-release N` on `annotate` and `identifiability`.

```python
server = ensembl.resolve_server(110)   # 'https://jul2023.rest.ensembl.org' -- GENCODE 44
cfg = annotate.build_config("LEPR", release=110)
res = identifiability.analyze(cfg, ensembl_release=110)
res["annotation"]                      # {'ensembl_release': 110, 'fetched_release': 110}
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

**`identifiability.kmers(seq, k) -> set`**
Return the set of length-`k` substrings (k-mers) of `seq` (upper-cased).

**`identifiability.analyze(config, k=31, sequences=None, *, canonical=True, window=None, background_sequences=None, background_fasta=None, background_gene_transcripts="auto", species=None, read_length=100, frag_mean=200.0, frag_sd=60.0, paired=True, depth=30_000_000, mean_efflen=1500.0, tpm=10.0, n_donors=1, conditioning_tau=10.0, min_informative_reads=50.0, min_log2fc=None, ensembl_release=None, inputs_out=None, retries=5, retry_wait=1.0) -> dict`**
The three layers the `identifiability` command reports -- sequence uniqueness, the read
model, and estimability on the class-collapsed compatibility system -- for each group and
for the contrast named in `primary_comparison`. `sequences` and `background_sequences`
are optional `{transcript_id: cdna}` dicts; whatever is missing is fetched from Ensembl
release `ensembl_release` (None: the current one), and nothing is fetched when everything
was supplied. Pass a dict as `inputs_out` to get back the sequence the run used
(`sequences`, `background_sequences`, `fetched_release`) and the system it was built at
(`k`, `window`, `canonical`, none of which is in the config), which `io.save_inputs`
writes to a file. Returns a dict with `k`, `window`, `canonical`; `annotation`
(`ensembl_release`, `fetched_release`); `background`; `design`; `groups` (per group:
`verdict`, `reasons`, `n_unique_kmers`, `unique_length`, `n_blocks`,
`expected_informative_reads`, `estimable`, `conditioning_factor`, `gls_relative_se`,
`min_resolvable_log2fc`, `coherence`, ...); `contrast` (the same estimability fields for
the class contrast, plus the effective-length and distinguishing-window summaries);
`counting_noise`; `gene_total`; `effect_resolvable` (None when no `min_log2fc` was
given); `n_compatibility_classes`; `verdict` and `reasons`. The docstring defines each.

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

A group with zero unique k-mers is not separable by short reads and is flagged
(`distinguishable: False`). Neither this flag nor `verdict` sets the CLI's exit status:
see the exit codes in the README.

---

## `extract`

**`extract.extract(config, quantdir, samplemap, cohort) -> list`**
Read `quantdir/<donor>/quant.sf` Salmon outputs and sum TPM per isoform group.
Returns `[(donor, condition, {group: tpm}), ...]`. Raises `FileNotFoundError` if no
`quant.sf` is found and `ValueError` if a file lacks the `Name`/`TPM` columns.

**`extract.run(config, quantdir, samplemap, cohort, out) -> int`**
Run `extract` and write a per-donor CSV to `out`; returns the number of donors.

---

## `stats`

**`stats.paired_stat(A, B) -> (n, n_A>B, P, median_fold)`**
Donor-level two-sided exact Wilcoxon signed-rank test (`scipy.stats.wilcoxon`) plus
median fold-change. Handles edge cases: empty input → NaNs; all-tied pairs → P is
NaN (test undefined); non-finite ratios are dropped from the fold-change.

**`stats.run(config, condition, cohorts, out) -> dict`**
`cohorts` is `{name: perdonor.csv}`. Writes `<out>.{png,pdf,svg}` and
`<out>_stats.csv`, and returns `{"per_cohort": [...], "combined": (n, n_gt, P, fold)}`.

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
**`io.save_inputs(path, captured, config, release, version, background_fasta=None) -> dict`** — write the sequence an `identifiability` run used (`captured`, from `analyze(..., inputs_out=...)`) with the release it came from, as `"format": "isoform-dominance/inputs/1"` (`io.INPUTS_FORMAT`). `analysis` records the `k`, `window` and `canonical` of the run, none of which the config holds. A `background_fasta` is recorded by path, size and SHA-256, not copied. Backs `--save-inputs`.
**`io.load_inputs(path) -> dict`** — read a file `save_inputs` wrote; `io.InputError` if it is not one, and for every field the caller goes on to read: `sequences` and `background_sequences` as objects of id to sequence, an `ensembl_release` that is a number or null, an `analysis` with integer `k` and `window` and boolean `canonical`, and a `background_fasta` that is null or has a path and a sha256. Which transcripts are a class and which are background is decided by the config the rerun is given, not by the saved grouping. Backs `--inputs`.
**`io.file_sha256(path) -> str`** — hex SHA-256 of a file's bytes.
**`io.transcript_to_group(groups) -> dict`** — invert `{group: [ENST...]}` to `{ENST(no version): group}`.
**`io.load_sample_map(path) -> dict`** — read a `donor,condition[,SRR]` CSV to `{donor: condition}`.
**`io.primary_pair(config) -> (gA, gB)`** — the two groups named in `primary_comparison`.
