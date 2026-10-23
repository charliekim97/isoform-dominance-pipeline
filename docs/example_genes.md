# Example genes

LEPR, NTRK2, NTRK3 and FLT1, each at GENCODE 44 (Ensembl 110) and GENCODE 50 (Ensembl 116),
from the configs `annotate --gtf` proposed and the sequence `identifiability --save-inputs`
saved, both kept in `example/saved/` (its README says how they were made). Run the commands
in a clone of the repository, from its root: `example/saved/` is not part of the installed
package, and none of them makes a network request. The `annotate` command below reads
GENCODE 50's GTF as well, which is not in the repository.

Each row is `isoform-dominance identifiability --config example/saved/<gene>.gencode<N>.config.json
--inputs example/saved/<gene>.gencode<N>.inputs.json` at the default design: the proposed
comparison with each class's transcript count, the `VERDICT`, and the `min |log2FC|` and
`conditioning` of the contrast line; and the exit status of the same command with
`--min-log2fc 0.5`. These are the fragment estimand's, which decides by default. The last two
columns are the abundance estimand's, the length-corrected class abundance, whose class ratio
is the ratio that summed TPM estimates (README, "Two estimands"): the contrast's
`min |log2FC|`, and the exit status with `--min-log2fc 0.5 --estimand abundance`.
`tests/test_docs_examples.py` recomputes every cell.

| Gene | GENCODE (Ensembl) | Proposed comparison (transcripts) | Verdict | Contrast min \|log2FC\| | Contrast conditioning | Exit status, `--min-log2fc 0.5` | Abundance: contrast min \|log2FC\| | Abundance: exit status |
|---|---|---|---|---|---|---|---|---|
| LEPR | 44 (110) | `iso_906aa` (2) vs `iso_1165aa` (1) | `weakly_identifiable` | 0.318 | 12.82 | 0 | 0.312 | 0 |
| LEPR | 50 (116) | `iso_896aa` (7) vs `iso_1165aa` (2) | `weakly_identifiable` | 0.234 | 117.68 | 0 | 0.225 | 0 |
| NTRK2 | 44 (110) | `iso_477aa` (18) vs `iso_822aa` (10) | `not_identifiable` | n/a | inf | 3 | 0.134 | 0 |
| NTRK2 | 50 (116) | `iso_477aa` (22) vs `iso_822aa` (12) | `not_identifiable` | n/a | inf | 3 | 0.133 | 0 |
| NTRK3 | 44 (110) | `iso_617aa` (2) vs `iso_839aa` (6) | `identifiable` | 0.295 | 9.14 | 0 | 0.261 | 0 |
| NTRK3 | 50 (116) | `iso_604aa` (7) vs `iso_825aa` (20) | `not_identifiable` | n/a | inf | 3 | 0.068 | 0 |
| FLT1 | 44 (110) | `iso_733aa` (1) vs `iso_1338aa` (1) | `identifiable` | 0.272 | 4.64 | 0 | 0.272 | 0 |
| FLT1 | 50 (116) | `iso_733aa` (1) vs `iso_1338aa` (6) | `weakly_identifiable` | 0.327 | 181.40 | 0 | 0.328 | 0 |

The verdict is a screening flag, and the exit status answers `--min-log2fc`; the
[NTRK2 walkthrough](tutorial_NTRK2.md) runs NTRK2 at GENCODE 50 and quotes what the command
says of both.

## NTRK3

The proposed comparison is not the same at the two releases: `iso_617aa` against
`iso_839aa` at GENCODE 44, `iso_604aa` against `iso_825aa` at GENCODE 50. The first is
`identifiable`, the second `not_identifiable`, and the exit status at `--min-log2fc 0.5`
goes from 0 to 3. As the abundance estimand (the last two columns) the contrast is estimable at
both releases, and the exit status is 0 at both: at GENCODE 50 no direction the system cannot
see changes either class's abundance, while some change a class's fragment total, by moving
abundance between transcripts of the class that differ in length. At GENCODE 44 `iso_617aa`
was itself chosen by a tie (below).

## FLT1, and proposals decided by a tie

At both releases FLT1's alternative class is one transcript, and so are two other clusters,
and `annotate` says so:

```bash
isoform-dominance annotate --gene FLT1 --gtf gencode.v50.annotation.gtf.gz \
    --out FLT1.gencode50.config.json
# Proposed groups for FLT1 -> FLT1.gencode50.config.json
#   iso_1338aa: 6 transcripts
#   iso_733aa: 1 transcript
#   primary_comparison: ['iso_733aa', 'iso_1338aa']
#   REVIEW _proposed/_clusters and rename groups before use.
# exit status 0
```

The config it wrote is kept as `example/saved/FLT1.gencode50.config.json`. On standard error:

> NOTE: iso_733aa was chosen by a tie. 2 other cluster(s) have as many transcripts (1):
> 687 aa at acceptor 28390104, 665 aa at acceptor 28389632. The rule is the non-canonical
> cluster with the most transcripts; a tie goes to the longer representative protein, and a
> tie on both to the lower terminal-acceptor coordinate, so a release that adds one
> transcript to either can change the proposed comparison. Pick the one you mean; each is
> listed under _clusters.

At GENCODE 50 the contrast's conditioning factor, 181.40, is above the default `--tau` of
10, which makes the verdict `weakly_identifiable`; at GENCODE 44 it is 4.64. Each config
records the clusters a tie was decided against under `_proposal.tied_with`:

```python
from isoform_dominance import io

for name in ("NTRK3.gencode44", "FLT1.gencode44", "FLT1.gencode50", "LEPR.gencode44"):
    cfg = io.load_config("example/saved/%s.config.json" % name)
    alt = cfg["primary_comparison"][0]
    tied = [(t["rep_protein_aa"], t["transcripts"]) for t in cfg["_proposal"]["tied_with"]]
    print(name, alt, len(cfg["groups"][alt]), "tied with", tied)
print(io.load_config("example/saved/LEPR.gencode44.study.config.json")["groups"])
# NTRK3.gencode44 iso_617aa 2 tied with [(612, ['ENST00000317501', 'ENST00000695463'])]
# FLT1.gencode44 iso_733aa 1 tied with [(687, ['ENST00000615840']), (665, ['ENST00000639477'])]
# FLT1.gencode50 iso_733aa 1 tied with [(687, ['ENST00000615840']), (665, ['ENST00000639477'])]
# LEPR.gencode44 iso_906aa 2 tied with [(896, ['ENST00000371060', 'ENST00000616738'])]
# {'short': ['ENST00000371060', 'ENST00000616738'], 'long': ['ENST00000349533']}
```

## LEPR

At GENCODE 44 the 906 aa and 896 aa clusters have two transcripts each, and the longer
protein is proposed; the README's first item ("What makes it more than a quantifier") gives
the release history. The cluster it was chosen over, ENST00000371060 and ENST00000616738, is
the short class of the choroid-plexus study, `short` in
`example/saved/LEPR.gencode44.study.config.json` (last line above); `isoform-dominance
selftest` runs those classes on `LEPR.gencode44.inputs.json`, from copies the package
carries. At GENCODE 50 the 896 aa class has seven transcripts and is the one proposed.
