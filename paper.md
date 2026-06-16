---
title: 'isoform-dominance: a Python command-line tool for functional isoform-class dominance from bulk RNA-seq, with a short-read identifiability guardrail'
tags:
  - Python
  - RNA-seq
  - alternative splicing
  - transcript usage
  - isoform quantification
  - reproducibility
authors:
  - name: Sangeon Kim
    orcid: 0009-0005-7227-323X
    affiliation: 1
affiliations:
  - name: Division of Endocrinology, Diabetes and Metabolism, Beth Israel Deaconess Medical Center, Harvard Medical School, Boston, MA, USA
    index: 1
date: 14 June 2026
bibliography: paper.bib
---

# Summary

Many biological questions reduce to a single, focused comparison: *for one gene,
which functional isoform class predominates in a given tissue or condition?* A
common example is a receptor with a signalling-competent long isoform and a
truncated short isoform that retains ligand binding but loses the intracellular
effector domain. Answering this from bulk RNA-seq requires four steps that are
usually stitched together by hand: deciding which annotated transcripts belong to
each functional class, checking whether short reads can even distinguish those
classes, quantifying each class per sample, and applying a paired test across
donors.

`isoform-dominance` packages these four steps behind a single Python
command-line interface. Given a gene symbol it proposes isoform groups from
Ensembl by clustering protein-coding transcripts on their 3' terminal-exon splice
acceptor [@frankish2021]; it then reports, per group, the number of unique
*k*-mers in the group's cDNA and refuses to treat a group as measurable when it
has no unique sequence — an explicit short-read **identifiability** check.
Quantification is read from `Salmon` `quant.sf` files [@patro2017] and summarised
per donor; group dominance is tested with a donor-level two-sided exact Wilcoxon
signed-rank test [@virtanen2020], reported per cohort and combined, alongside a
median fold-change effect size. An optional contamination control tests whether a
dominance signal is driven by cell-type marker contamination rather than true
isoform usage. A bundled, download-free self-test reproduces a published result
(short-isoform dominance of the leptin receptor in human choroid plexus across
two independent cohorts) on a clean machine in seconds.

# Statement of need

Researchers studying a specific gene often need to answer one focused question:
*which functional isoform class predominates in this tissue or condition?* This
arises whenever a gene produces functionally divergent isoforms through an
alternative terminal exon — for example a signalling-competent long receptor and a
truncated or secreted short form. The target audience is bench and computational
biologists testing such a hypothesis from bulk RNA-seq, who do not need a
genome-wide screen but do need the focused comparison done correctly and
reproducibly.

Doing this by hand has a recurring, often-unstated pitfall: short-read
quantifiers can only apportion reads to an isoform class that carries *unique*
sequence, so a reported fold-change between two classes can be an artefact when the
classes are not actually separable by short reads. `isoform-dominance` addresses
this by going directly from a *gene symbol* to a reviewed isoform-group proposal,
making the separability assumption explicit and machine-checkable *before* any
quantification, and packaging the comparison as a scriptable CLI with a
reproducible, download-free self-test suitable for inclusion alongside a
manuscript.

# State of the field

Differential transcript usage (DTU) is a mature area with strong tools, but they
target a different task and audience. `DEXSeq` [@anders2012],
`DRIMSeq` [@nowicka2016] and `satuRn` [@gilis2022] perform genome-wide DTU testing
and assume the user has already produced a transcript-by-sample count matrix and
defined the transcript groups. `IsoformSwitchAnalyzeR` [@vittingseerup2019] adds
rich functional annotation of isoform switches but is an R/Bioconductor workflow
in which transcript grouping and import are configured by the analyst.
`fishpond`/`swish` [@zhu2019] rigorously propagates quantification uncertainty
using Salmon [@patro2017] inferential replicates, but likewise operates on a
prepared dataset. Transcript-level import for these workflows is typically handled
by `tximport` [@soneson2016].

**Build vs. contribute.** These tools solve genome-wide discovery and were not the
right place to add this functionality: none takes a bare gene symbol as input, and
none performs an up-front, short-read *identifiability* check on the proposed
grouping. `isoform-dominance` is deliberately narrow and complementary rather than
competing — it does not attempt genome-wide discovery, delegates quantification to
Salmon rather than re-implementing it, and presents group proposals for review
rather than as final calls. Its distinct scholarly contribution is the end-to-end,
reproducible packaging of the focused dominance question and, in particular, the
identifiability guardrail, which prevents the silent failure mode of reporting a
fold-change between isoform classes that short reads cannot separate. Where prior
similar functionality exists (e.g. uncertainty propagation in `fishpond`), it is
cited rather than reinvented.

The tool generalises across unrelated biology. Beyond the bundled leptin-receptor
example, `annotate` + `identifiability` cleanly recover, and confirm as short-read
separable: the full-length vs truncated kinase receptors of *NTRK2* (TrkB-FL vs
TrkB.T1) and *NTRK3* (TrkC), and the membrane vs secreted-decoy forms of *FLT1*
(full VEGFR1 vs soluble sFlt-1) — worked examples included in the documentation.

# Software design

The package is organised as small, independently testable modules behind a single
`argparse` CLI: `annotate` (Ensembl REST query and terminal-exon clustering),
`identifiability` (group-unique *k*-mer counting), `extract` (Salmon `quant.sf`
parsing and per-donor group TPM), `stats` (paired exact Wilcoxon and plotting),
and `contamination` (marker-correlation control). Configuration is a single JSON
file describing the isoform groups, so the same pipeline applies to any gene
without code changes. Network-dependent steps (`annotate`, `identifiability`)
accept offline inputs to keep the test suite deterministic and to allow use on
compute nodes without internet access. The design deliberately delegates the
heavy, well-solved problem of transcript quantification to Salmon rather than
re-implementing it, and focuses original effort on the grouping, identifiability,
and reproducible-comparison layers. The annotation and sequence queries use the
Ensembl REST API [@yates2015]; numerical work uses NumPy [@harris2020] and SciPy
[@virtanen2020], and figures use Matplotlib [@hunter2007]. Dependencies are kept
minimal so the package is straightforward to install and audit.

# Research impact

`isoform-dominance` is the generalized, reusable implementation of an
isoform-usage analysis developed for a study of leptin-receptor isoforms in the
human choroid plexus. Rather than aspirational claims, concrete and reproducible
evidence is provided. First, the bundled, download-free self-test reproduces that
study's published result — short-isoform dominance across two independent public
cohorts [@gse228458; @gse137619]; combined n = 11, paired Wilcoxon
P = 9.8×10⁻⁴ — on a clean machine in seconds, and the package is archived on
Zenodo with a citable DOI. Second, the study's standalone reproducibility
repository runs its analysis through this package, so the published result is
regenerated by the software itself rather than by separate one-off code. Third,
the method's generality is shown on three further, biologically unrelated loci
with verifiable outputs (the *NTRK2*, *NTRK3*, and *FLT1* worked examples in the
documentation), each confirming the proposed isoform classes are short-read
separable. *The citation for the associated leptin-receptor publication will be
added here once it is available.* The tool targets researchers asking which
functional isoform class predominates in a tissue or condition.

# AI usage disclosure

Generative AI assistance (Anthropic Claude) was used during development for code
scaffolding, refactoring, test drafting, documentation, and copy-editing of this
paper. All AI-assisted outputs were reviewed, edited, and validated by the author,
who made the core design decisions (the terminal-exon grouping approach, the
identifiability guardrail, the statistical model, and the package architecture)
and is fully responsible for the correctness, originality, and licensing of the
submitted materials.

# Acknowledgements

Computation used the O2 High Performance Compute Cluster at Harvard Medical
School. The author thanks the leptin-receptor study team for the data that
motivated the bundled example.

# References
