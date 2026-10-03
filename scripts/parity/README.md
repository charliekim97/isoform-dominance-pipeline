# Parity of the local-annotation mode with REST

`annotate --gtf` and `identifiability --gtf --transcripts-fasta` (issue #13) are meant to give
what `annotate` and `identifiability` give from Ensembl REST at the release the files are
of. The tests hold that on a seven-gene extract of GENCODE 50 (`tests/data/gencode_mini`), and
a weekly job holds the extract to live REST 116 (`.github/workflows/ensembl-nightly.yml`,
`gtf-parity`). These scripts make the full comparison, on whole files, by hand. They take
every path as an argument and are not run in CI.

| script | what it does | network |
|---|---|---|
| `record_rest.py` | records what REST answers for a list of genes at one release: `lookup/symbol`, `xrefs/symbol`, `lookup/id?expand=1` of every gene they name, and the MD5 of every transcript's cDNA, in the format of `tests/data/gencode_mini/rest116_mini.json.gz` | yes |
| `compare.py` | compares, gene by gene, the config `annotate` writes from those answers with the one `annotate --gtf` writes from the GTF (`annotation_source` aside), and every transcript of the gene the GTF chose, from the FASTA, with REST's cDNA by MD5; exit 0 when all are equal | no |
| `timing.py` | times `annotate --gtf` (exact name, then the name in lower case) and `identifiability --gtf --transcripts-fasta` through the command line | no |

## The full comparison

GENCODE 44, 48 and 50 are Ensembl 110, 114 and 116. With `genes.txt` the 109-gene survey
panel, one symbol per line:

```sh
for pair in 44:110 48:114 50:116; do
  g=${pair%%:*}; e=${pair##*:}
  python scripts/parity/record_rest.py --release $e --genes genes.txt --out rest$e.json.gz
  python scripts/parity/compare.py --gtf gencode.v$g.annotation.gtf.gz \
      --fasta gencode.v$g.transcripts.fa.gz --rest rest$e.json.gz --genes genes.txt \
      --json parity$e.json
done
python scripts/parity/timing.py --gtf gencode.v50.annotation.gtf.gz \
    --fasta gencode.v50.transcripts.fa.gz --gene LEPR
```

`record_rest.py` needs the release to be served: by `rest.ensembl.org` for 116, by the REST
archive for 110 and 114 while Ensembl keeps them. `compare.py` first cuts the lines of the
genes asked about from the GTF, and their records from the FASTA, in one pass each (as the
extract in `tests/data/gencode_mini` was cut); `--whole` reads the whole files for every
gene instead.

What the design study measured, with the prototype of this reader: every one of the 109
configs equal at each release, and the cDNA byte-identical for 2005, 2066 and 4396
transcripts at 110, 114 and 116. `compare.py` counts every transcript of each gene the GTF
chose that REST's answers hold, so its count is that of the genes recorded.

On the extract (`tests/test_parity_scripts.py` runs this):

```sh
python scripts/parity/compare.py --gtf tests/data/gencode_mini/gencode.v50.mini.gtf.gz \
    --fasta tests/data/gencode_mini/gencode.v50.mini.transcripts.fa.gz \
    --rest tests/data/gencode_mini/rest116_mini.json.gz
# release 116: configs equal 6/6 (annotation_source aside); cDNA byte-identical 217/217; ...
```

HERC3, the seventh gene of the extract, has no recorded REST answer: its two genes stop
the GTF path with `AmbiguousGene` (`tests/test_gtf_parity.py`), as two such genes stop REST's
(`tests/test_gene_choice.py`).
