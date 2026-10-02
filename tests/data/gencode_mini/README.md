# GENCODE 50 and Ensembl 116 mini fixtures (local annotation, issue #13)

For the tests of the local-annotation mode. They were cut on 2026-09-29 from the full release files and from REST
responses fetched 2026-09-23 to 2026-09-29.

| file | contents | SHA-256 |
|---|---|---|
| `gencode.v50.mini.gtf.gz` | the header and every line of GENCODE 50's comprehensive `gencode.v50.annotation.gtf.gz` (dated 2026-04-08) for the genes named LEPR, FOXO1, STK11, AXIN1, GSK3B, CD99 (chrX ENSG00000002586, chrY ENSG00000292348) and HERC3 (ENSG00000138641, ENSG00000287542): 6,334 lines | `0797a95dbf00ccd0a511804937ef632cc03e09cd4a898eb603ab29fe40d76e60` |
| `gencode.v50.mini.transcripts.fa.gz` | the records of `gencode.v50.transcripts.fa.gz` for those genes' transcripts, without the `_PAR_Y` records: 306 | `83a71068b0844db4eab3b17d2e4d2418ba15db1839df344fbcdb41b40c2ac91f` |
| `rest116_mini.json.gz` | Ensembl 116 REST (`rest.ensembl.org`): `lookup/id/<gene>?expand=1` for LEPR, FOXO1, STK11, AXIN1, GSK3B and both CD99 genes, `lookup/symbol` for CD99, and the MD5 of 217 transcripts' cDNA, trimmed to the fields the package reads | `30a0d430e15ccace6f133267deff04df3471117db1480082199a34da0957970c` |

GENCODE 50 is Ensembl 116.
