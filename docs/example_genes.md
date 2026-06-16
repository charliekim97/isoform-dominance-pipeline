# Example genes: generality across unrelated biology

The bundled self-test uses the leptin receptor (LEPR), and
[`tutorial_NTRK2.md`](tutorial_NTRK2.md) walks through NTRK2 step by step. To show
that the `annotate` → `identifiability` workflow is not tuned to one gene or one
kind of biology, here are the verified outputs for several functionally split
genes spanning three different mechanisms. All values come from running the tool
live against the Ensembl REST API; transcript and *k*-mer counts reflect the
annotation at the time of writing and will drift as the annotation updates, but the
two-class structure and short-read separability are stable.

| Gene | Functional split (mechanism) | Long / full class | Short / alt class | Both short-read separable? |
|---|---|---|---|---|
| **LEPR** | transport vs signalling (alt terminal exon) | LepRb (long, 1165 aa) | LepRa (short) | yes (bundled self-test) |
| **NTRK2** | kinase-competent vs truncated receptor | TrkB-FL (822 aa) | TrkB.T1 (477 aa) | yes — 6793 / 5639 unique 31-mers |
| **NTRK3** | kinase-competent vs truncated receptor | TrkC (825 aa) | truncated (604 aa) | yes — 18442 / 2157 unique 31-mers |
| **FLT1** | membrane receptor vs secreted decoy | VEGFR1 (1338 aa) | sFlt-1 (soluble, 733 aa) | yes — 5656 / 568 unique 31-mers |

These span a neurotrophin-kinase truncation (NTRK2/NTRK3), an angiogenic
soluble-decoy receptor produced by an alternative terminal exon (FLT1 / sFlt-1,
central to pre-eclampsia biology), and the bundled leptin-receptor transport case
— biology with nothing in common beyond the structural motif the tool targets.

## Reproducing the FLT1 example

```bash
isoform-dominance annotate --gene FLT1 --out flt1.json
#   iso_1338aa: 6 transcripts   (full membrane VEGFR1)
#   iso_733aa : 1 transcript    (soluble sFlt-1 decoy)

isoform-dominance identifiability --config flt1.json
#   [OK] iso_1338aa: 5656 unique k-mers
#   [OK] iso_733aa :  568 unique k-mers
#   primary_comparison distinguishable by short reads: True
```

Note that the soluble class (sFlt-1) carries far fewer unique *k*-mers (568) than
the full receptor — it shares most of its sequence with the membrane form and
differs mainly in its alternative terminal exon. The guardrail still passes here
(568 > 0), but this is exactly the kind of case where the identifiability check
earns its keep: a class that is *barely* separable is flagged quantitatively
rather than assumed.

After reviewing and renaming the groups in the JSON (e.g. `VEGFR1` / `sFlt1`), the
remaining `extract` → `stats` (→ `qc`) steps are identical to the LEPR workflow.
