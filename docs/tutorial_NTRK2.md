# Worked example: NTRK2 (a gene unrelated to the bundled data)

The bundled self-test uses the leptin receptor (LEPR). To show that the workflow
is general and not tuned to one gene, here is the same `annotate` →
`identifiability` flow applied to **NTRK2**, which encodes the neurotrophin
receptor TrkB. NTRK2 is a textbook case of a functionally split gene: a
full-length, kinase-competent receptor (**TrkB-FL**) and a truncated isoform
class (**TrkB.T1**) that keeps the extracellular ligand-binding region but lacks
the intracellular tyrosine-kinase domain.

This example requires network access (it queries the Ensembl REST API).

## 1. Propose isoform groups from the gene symbol

```bash
isoform-dominance annotate --gene NTRK2 --out ntrk2.json
```

Observed output:

```
Proposed groups for NTRK2 -> ntrk2.json
  iso_822aa: 12 transcripts
  iso_477aa: 22 transcripts
  primary_comparison: ['iso_477aa', 'iso_822aa']
  REVIEW _proposed/_clusters and rename groups before use.
```

The two clusters correspond to the long kinase-competent receptor (~822 aa,
TrkB-FL) and the truncated class (~477 aa, TrkB.T1). As always, the proposal is
for review: open `ntrk2.json` and rename `iso_822aa` / `iso_477aa` to
`TrkB_FL` / `TrkB_T1` before downstream use, and confirm the transcript
membership matches your intended definition.

## 2. Check short-read identifiability *before* quantifying

```bash
isoform-dominance identifiability --config ntrk2.json
```

Observed output:

```
Identifiability (k=31):
  [OK] iso_822aa: 6793 unique k-mers (12 transcripts)
  [OK] iso_477aa: 5639 unique k-mers (22 transcripts)
  primary_comparison distinguishable by short reads: True
```

Both groups carry thousands of unique 31-mers, so short-read quantification can
separate TrkB-FL from TrkB.T1. Had one group contained no unique sequence, the
command would have flagged it (`NOT DISTINGUISHABLE`) and exited non-zero,
stopping you from reporting a fold-change the data cannot support.

## 3. From here

With a reviewed `ntrk2.json` and Salmon `quant.sf` files for your samples, the
remaining steps are identical to the LEPR workflow in the
[README](../README.md): `extract` → `stats` (→ optional `qc`).

> Transcript counts and *k*-mer totals above reflect the Ensembl annotation at
> the time of writing and will change as the annotation is updated; the
> interpretation (two clearly distinguishable functional classes) is stable.
