# `identifiability --json` of the 2.7.0 code on the saved example inputs

`tests/test_saved_identifiability.py` runs the current code on the same inputs and compares its report with these
files: the same keys, the same strings, booleans, integers and nulls, and every real number within a relative
tolerance of 1e-9 (absolute 1e-12), after removing the one key 2.8.0 adds to the report, `unique_read_counting_se`.
The tolerance is there because a report holds about 1e-16 `residual` values whose last digits move with the BLAS
build; nothing else about the report may move.

They were written on 2026-10-17 by the `src/` of the 2.7.0 release (commit 4350016; the `src/` of the commit that
opened 2.8.0.dev0 differs from it only in the version string), from the repository root, with the saved inputs of
`example/saved/` (its README says how they were made), the default design and no `--min-log2fc`:

```bash
for G in FLT1 LEPR NTRK2 NTRK3; do for N in 44 50; do
  isoform-dominance identifiability --config example/saved/$G.gencode$N.config.json \
      --inputs example/saved/$G.gencode$N.inputs.json --json > $G.gencode$N.json
done; done
isoform-dominance identifiability --config example/saved/LEPR.gencode44.study.config.json \
    --inputs example/saved/LEPR.gencode44.inputs.json --json > LEPR.gencode44.study.json
```

Every run exited 0 and wrote nothing to stderr. The report does not depend on the paths given: the same files came
out of a run from another directory with absolute paths.

| file | `verdict` | contrast estimable | contrast `min_resolvable_log2fc` |
|---|---|---|---|
| `FLT1.gencode44.json` | identifiable | yes | 0.2717 |
| `FLT1.gencode50.json` | weakly_identifiable | yes | 0.3270 |
| `LEPR.gencode44.json` | weakly_identifiable | yes | 0.3183 |
| `LEPR.gencode50.json` | weakly_identifiable | yes | 0.2342 |
| `NTRK2.gencode44.json` | not_identifiable | no | null |
| `NTRK2.gencode50.json` | not_identifiable | no | null |
| `NTRK3.gencode44.json` | identifiable | yes | 0.2949 |
| `NTRK3.gencode50.json` | not_identifiable | no | null |
| `LEPR.gencode44.study.json` | identifiable | yes | 0.1290 |
