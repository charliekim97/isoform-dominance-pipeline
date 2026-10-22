# Worked example: NTRK2

The bundled self-test uses the leptin receptor (LEPR). This page runs the same `annotate` →
`identifiability` flow on NTRK2 and shows what each command prints. Run the commands in a
clone of the repository, from its root: they read the inputs saved in `example/saved/`,
which the installed package does not include, and none of them makes a network request.

## 1. The proposed classes (GENCODE 50 = Ensembl 116)

`annotate --gtf` proposes the classes from GENCODE's comprehensive GTF:

```bash
isoform-dominance annotate --gene NTRK2 --gtf gencode.v50.annotation.gtf.gz \
    --out NTRK2.gencode50.config.json
# Proposed groups for NTRK2 -> NTRK2.gencode50.config.json
#   iso_822aa: 12 transcripts
#   iso_477aa: 22 transcripts
#   primary_comparison: ['iso_477aa', 'iso_822aa']
#   REVIEW _proposed/_clusters and rename groups before use.
# exit status 0
```

The GTF is not in the repository. The config this command wrote is
`example/saved/NTRK2.gencode50.config.json`, and `example/saved/README.md` says how it and the
saved sequence were made. A class is the protein-coding transcripts that share a 3'
terminal-exon splice acceptor, named by the length of its representative protein; the
config lists each class's transcript ids, and the other clusters under `_clusters`. Review
and rename the classes before use.

## 2. Identifiability on the saved sequence

```bash
isoform-dominance identifiability --config example/saved/NTRK2.gencode50.config.json \
    --inputs example/saved/NTRK2.gencode50.inputs.json
# Identifiability (window=31, canonical k-mers)
#   annotation: Ensembl release 116; sequence from saved inputs (release 116, example/saved/NTRK2.gencode50.inputs.json)
#   background: 30 same-gene transcript(s)
#   design: paired 100bp reads, fragments 200+-60, depth 30M, TPM 10, n=1
#   [not_identifiable] iso_822aa: min |log2FC| n/a; class 31 unique k-mers; best transcript 57 bp in 1 block(s), ~34 informative reads; conditioning inf
#   [not_identifiable] iso_477aa: min |log2FC| n/a; class 287 unique k-mers; best transcript 60 bp in 1 block(s), ~35 informative reads; conditioning inf
#   contrast iso_477aa vs iso_822aa: min |log2FC| n/a; not_identifiable, conditioning inf
#   effective length: iso_477aa 6232 bp vs iso_822aa 7940 bp (mean per transcript), log2 ratio -0.35
#   distinguishing windows (0 = 5' end, 1 = 3' end of each transcript): iso_477aa median 0.24 [0.05-0.99], 360 positions over 22 transcript(s); iso_822aa median 0.04 [0.04-1.00], 39 positions over 12 transcript(s)
#   unique-read counting SE (best transcript per class) on log2 ratio: 0.347 per donor, min resolvable |log2FC| 0.680 at n=1
#   VERDICT: not_identifiable
# exit status 0
```

Without `--min-log2fc` the exit status is 0 whatever the verdict. Given the smallest effect
you need, the command answers that question in its exit status:

```bash
isoform-dominance identifiability --config example/saved/NTRK2.gencode50.config.json \
    --inputs example/saved/NTRK2.gencode50.inputs.json --min-log2fc 0.5
# Identifiability (window=31, canonical k-mers)
#   annotation: Ensembl release 116; sequence from saved inputs (release 116, example/saved/NTRK2.gencode50.inputs.json)
#   background: 30 same-gene transcript(s)
#   design: paired 100bp reads, fragments 200+-60, depth 30M, TPM 10, n=1
#   [not_identifiable] iso_822aa: min |log2FC| n/a; class 31 unique k-mers; best transcript 57 bp in 1 block(s), ~34 informative reads; conditioning inf
#   [not_identifiable] iso_477aa: min |log2FC| n/a; class 287 unique k-mers; best transcript 60 bp in 1 block(s), ~35 informative reads; conditioning inf
#   contrast iso_477aa vs iso_822aa: min |log2FC| n/a; not_identifiable, conditioning inf
#   effective length: iso_477aa 6232 bp vs iso_822aa 7940 bp (mean per transcript), log2 ratio -0.35
#   distinguishing windows (0 = 5' end, 1 = 3' end of each transcript): iso_477aa median 0.24 [0.05-0.99], 360 positions over 22 transcript(s); iso_822aa median 0.04 [0.04-1.00], 39 positions over 12 transcript(s)
#   unique-read counting SE (best transcript per class) on log2 ratio: 0.347 per donor, min resolvable |log2FC| 0.680 at n=1
#   EFFECT SIZE: |log2FC| 0.500 NOT resolved at this design (iso_477aa, iso_822aa, the contrast beyond the linearisation limit, relative SE > 0.3)
#   VERDICT: not_identifiable
# exit status 3
```

## 3. What the command says of its verdict

The blocks above are standard output. On standard error the first run also says, word for
word:

> The verdict is a screening flag, not the exit status. It says whether each class total
> and the contrast lie in the row space of a sequence-derived compatibility surrogate at
> window=31, for exactly the transcripts in this config. It does not say whether a
> quantifier will recover the comparison, and it changes with the annotation release those
> transcripts came from. Pass --min-log2fc <the smallest effect you need> for an exit
> status on the resolvable effect size.

> Outside the row space of the compatibility surrogate at this window length: iso_822aa,
> iso_477aa, the class contrast. That surrogate is built from sequence, not from the
> observation model of a sequencing run, so this is a screening verdict rather than a
> statement about the data: a longer --window, a different grouping, or long reads may
> change it.

## 4. Unique k-mers between the two classes alone

The background is the gene's 30 other transcripts, which the inputs hold. Without it,
uniqueness is judged between the two classes only:

```bash
isoform-dominance identifiability --config example/saved/NTRK2.gencode50.config.json \
    --inputs example/saved/NTRK2.gencode50.inputs.json --no-gene-background
# Identifiability (window=31, canonical k-mers)
#   annotation: Ensembl release 116; sequence from saved inputs (release 116, example/saved/NTRK2.gencode50.inputs.json)
#   background: 0 same-gene transcript(s)
#   design: paired 100bp reads, fragments 200+-60, depth 30M, TPM 10, n=1
#   [not_identifiable] iso_822aa: min |log2FC| n/a; class 6793 unique k-mers; best transcript 6772 bp in 2 block(s), ~1369 informative reads; conditioning inf
#   [not_identifiable] iso_477aa: min |log2FC| n/a; class 5639 unique k-mers; best transcript 5356 bp in 2 block(s), ~1086 informative reads; conditioning inf
#   contrast iso_477aa vs iso_822aa: min |log2FC| n/a; not_identifiable, conditioning inf
#   effective length: iso_477aa 6232 bp vs iso_822aa 7940 bp (mean per transcript), log2 ratio -0.35
#   distinguishing windows (0 = 5' end, 1 = 3' end of each transcript): iso_477aa median 0.66 [0.48-0.84], 92484 positions over 22 transcript(s); iso_822aa median 0.63 [0.44-0.81], 72229 positions over 12 transcript(s)
#   unique-read counting SE (best transcript per class) on log2 ratio: 0.059 per donor, min resolvable |log2FC| 0.115 at n=1
#   VERDICT: not_identifiable
# exit status 0
```

6793 and 5639 are the two counts the version of this page before 2.9.0 showed, under
`primary_comparison distinguishable by short reads: True`: unique k-mers counted between
the configured classes alone. The verdict with them is the one above, `not_identifiable`.

## 5. GENCODE 44, and from here

The same comparison at GENCODE 44 (Ensembl 110), and LEPR, NTRK3 and FLT1 at both releases,
are in the [example genes](example_genes.md) table. With a reviewed config and Salmon
`quant.sf` files, the remaining steps are `extract` → `stats` (→ `qc`), as in the
[full workflow](index.md#full-workflow-real-data).
