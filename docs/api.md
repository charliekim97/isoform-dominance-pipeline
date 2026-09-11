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

---

## `annotate`

Propose isoform groups for a gene from Ensembl by clustering protein-coding
transcripts on their 3' terminal-exon splice acceptor.

**`annotate.fetch_transcripts(gene, species="homo_sapiens") -> dict`**
Query Ensembl and return `{gene, species, strand, transcripts: [{id, protein_aa, terminal_acceptor, is_canonical}]}`. Raises `ValueError` if no protein-coding transcripts with a translation are found. Network access required.

**`annotate.cluster_by_terminal_exon(info) -> list`**
Group the transcripts from `fetch_transcripts` by terminal-exon acceptor coordinate; returns clusters sorted by size, each `{acceptor, rep_aa, n, canonical, ids}`.

**`annotate.propose_groups(info) -> (groups, primary, clusters)`**
Choose the canonical cluster and the largest alternative cluster, returning `groups` ({label: [ids]}, labels like `iso_896aa`), the `primary_comparison` list (alternative first), and all clusters.

**`annotate.build_config(gene, species="homo_sapiens") -> dict`**
Convenience wrapper returning a complete, reviewable config dict (including `_proposed` notes and `_clusters`).

**`annotate.run(gene, out, species="homo_sapiens") -> dict`**
As `build_config`, but also writes the config JSON to `out`. Backs the `annotate` CLI subcommand.

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

**Retries.** Every Ensembl request (`annotate`'s lookup included) goes through
`ensembl.request`. It is retried on HTTP 429/500/502/503/504, connection errors and read
timeouts, and not on any other HTTP status. `retry` is `retries=` (default
`ensembl.DEFAULT_RETRIES = 5`, counted after the first attempt), `retry_wait=` (default
`ensembl.DEFAULT_RETRY_WAIT = 1.0` s before the first retry, doubling each time; a 429
waits for its `Retry-After` instead) and `timeout=` (default `ensembl.DEFAULT_TIMEOUT =
30` s per attempt). `annotate.build_config`/`run` and `identifiability.analyze` take the
same `retries`/`retry_wait`; on the CLI they are `--retries` and `--retry-wait`.

**`identifiability.kmers(seq, k) -> set`**
Return the set of length-`k` substrings (k-mers) of `seq` (upper-cased).

**`identifiability.analyze(config, k=31, sequences=None) -> dict`**
For each group, count k-mers unique to that group (not present in any other group).
`sequences` is an optional `{transcript_id: cdna}` dict; if omitted, sequences are
fetched from Ensembl. Returns:

```python
{
  "k": 31,
  "groups": {"<group>": {"n_unique_kmers": int, "n_transcripts": int, "distinguishable": bool}, ...},
  "primary_comparison": [...],
  "primary_distinguishable": bool   # True only if every primary group has >0 unique k-mers
}
```

A group with zero unique k-mers is not separable by short reads and is flagged
(`distinguishable: False`); the CLI exits non-zero in that case.

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

**`io.load_config(path) -> dict`** — load a config JSON.
**`io.transcript_to_group(groups) -> dict`** — invert `{group: [ENST...]}` to `{ENST(no version): group}`.
**`io.load_sample_map(path) -> dict`** — read a `donor,condition[,SRR]` CSV to `{donor: condition}`.
**`io.primary_pair(config) -> (gA, gB)`** — the two groups named in `primary_comparison`.
