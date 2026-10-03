# Outputs of the 2.6.0 code, for the invariants of 2.7.0

`tests/test_pinned_outputs.py` compares the current code's output with these files byte for byte.
They were written on 2026-10-03 at commit 6453060, whose `annotate.py`, `annotation_files.py`,
`extract.py` and `index_scope.py` are those of the 2.6.0 release (691962e); only the version
string differs.

| files | written by |
|---|---|
| `annotate_rest_<GENE>.json` | `annotate.run(gene, out)`, with the REST calls answered from `../gencode_mini/rest116_mini.json.gz` as `tests/test_gtf_parity.py` answers them |
| `annotate_gtf_<GENE>.json` | `annotate.run(gene, out, gtf="tests/data/gencode_mini/gencode.v50.mini.gtf.gz")` |
| `extract_selftest_<COHORT>.csv`, `.csv.index.json` | `extract.run(_selftest.CONFIG, ...)` on the self-test's synthetic quant.sf files (`_selftest.generate`) |

The genes are LEPR, FOXO1, STK11, AXIN1, GSK3B and CD99 (the chrX gene, by the pseudoautosomal
rule). The JSON files end without a newline, as `json.dump` leaves them.
