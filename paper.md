---
title: 'isoform-dominance: deciding whether a functional isoform-class comparison is measurable before running it'
tags:
  - Python
  - RNA-seq
  - alternative splicing
  - transcript usage
  - identifiability
  - reproducibility
authors:
  - name: Sangeon Kim
    orcid: 0009-0005-7227-323X
    affiliation: 1
affiliations:
  - name: Division of Endocrinology, Diabetes and Metabolism, Beth Israel Deaconess Medical Center, Harvard Medical School, Boston, MA, USA
    index: 1
date: 9 September 2026
bibliography: paper.bib
---

# Summary

Many biological questions reduce to one focused comparison: *for a given gene, which
functional isoform class predominates in a tissue or condition?* A common example is a
receptor's signalling-competent long isoform against a truncated short one that keeps ligand
binding but loses the intracellular effector domain. Answering it from bulk RNA-seq means
assigning transcripts to each class, quantifying each per sample and testing the contrast
across donors — and, first, establishing that short reads can resolve it.

`isoform-dominance` packages these steps behind one command-line interface. Given a gene
symbol it proposes isoform groups by clustering the annotation's protein-coding transcripts
[@frankish2021] on their 3' terminal-exon splice acceptor. Before quantification it answers
the measurability question in three layers: how much sequence is unique to each class
against the rest of the gene or a background FASTA; how many fragments would cover it at a
stated design, and hence the counting SE of the class ratio; and whether each class total,
and their contrast, is an *estimable* function of the fragment-class system. Quantification
is read from Salmon `quant.sf` files [@patro2017] per donor; the contrast is tested with a
donor-level Wilcoxon signed-rank test [@wilcoxon1945], per cohort with a bootstrap interval
and the finest p-value the design could resolve, and combined across cohorts. A bundled
self-test reproduces the reference result offline in seconds.

# Statement of need

The audience is biologists testing a hypothesis about one gene from bulk RNA-seq, who need
the focused comparison done correctly and reproducibly rather than a genome-wide screen.

Doing it by hand has an often-unstated pitfall. A short-read quantifier apportions fragments
among transcripts by solving a linear inverse problem, and a transcript's abundance — or a
group's — is recoverable only when it is an estimable function of that system. An EM
quantifier returns abundances for a non-identifiable model as readily as for an identifiable
one, stopping somewhere on the ridge where the likelihood is flat. A fold-change between two
classes the data cannot separate is then an artefact of tie-breaking, indistinguishable in
the output from a real one.

The obvious proxy — "does each class own a unique k-mer?", used up to v2.1.1 — is wrong in
both directions. A class whose only unique sequence is the thirty k-mers spanning one splice
junction passes while yielding a handful of informative fragments; a class nested inside
another owns no unique k-mer, yet its total is estimable, because the containing class is
pinned by its own unique sequence. Both failure modes are tests here.

`isoform-dominance` replaces the proxy with the estimability condition itself, reported
beside the smallest log2 fold change the stated design resolves and the expected number of
informative fragments. It withholds a pass when the contrast lies outside the row space, and
says "weakly identifiable, and here is why" where the honest answer is neither yes nor no.

# State of the field

**Differential transcript usage.** DTU tools target a different task. `DEXSeq` [@anders2012]
tests differential exon usage; `DRIMSeq` [@nowicka2016], `satuRn` [@gilis2022] and
`IsoformSwitchAnalyzeR` [@vittingseerup2019] test transcript-level DTU genome wide;
`fishpond`/`swish` [@zhu2019] propagates Salmon's quantification uncertainty into testing.
All begin from a count matrix, typically imported with `tximport` [@soneson2016], with the
groups already defined.

**Identifiability.** @hiller2009 gave criteria under which isoform abundances are uniquely
determined by RNA-seq observations; @ferrerbonsoms2022 made them read- and
fragment-length-aware and applied them genome wide. @zheng2022 found that for 20-47% of
Human Body Map transcripts the error from nonidentifiability is large enough that "the
ranking of expression between the transcript and other isoforms from the same gene cannot be
determined" — the comparison this package is asked to certify. `terminus` [@sarkar2020]
groups transcripts *post hoc, from the data*, so analysis proceeds at the resolution the
data support.

**The gap.** Those address transcript-level identifiability, or the resolution the data
permit. Here the classes are fixed *a priori, on biological grounds*, and the estimand is a
contrast of class totals. Neither a tool that begins from a count matrix nor a
transcript-level criterion has a place for that question, asked before quantification, so it
is answered here rather than by contributing to one. A class total is itself two estimands,
the expected fragment count `NumReads` sums and the length-corrected abundance whose class
ratio summed TPM estimates; they part when an invisible direction moves abundance between
transcripts of one class that differ in length. Since version 2.9 both are reported. The
fragment-total estimand decides by default, for continuity with earlier releases; the
abundance-total estimand, which matches what `extract` and `stats` compare, decides when
asked. Its precision figures were checked against Salmon's summed TPM in a simulation of 49
genes, by a gate fixed before its run, and did not pass, so the fragment-total estimand
stays the default: of the 7 contrasts estimable as abundance only against each gene's own
transcripts, they held for 6 and not for SMN2's. Found after the grade: against the
reference-chromosome transcripts of the simulation's index, SMN2's contrast is not estimable
as abundance, the transcripts of its paralogue SMN1 becoming columns.

Two questions hide under one word. A contrast is *estimable* when it lies in the row space
of the compatibility system, and *measurable* when it is recoverable at finite depth.
Transcript-level identifiability is not necessary for the first — two identical transcripts
are individually unidentifiable while their sum is estimable — nor does it discriminate when
it holds: full column rank makes *every* contrast estimable. It is silent on the second.
This package reports a verdict and a resolvable effect size together, on the user's own
contrast.

# Software design

Independently testable modules (`annotate`, `identifiability`, `extract`, `stats`,
`contamination`) sit behind one CLI; one JSON file names the groups, so any gene needs no
code changes, and every subcommand accepts `--json`.

*Uniqueness is judged against a background, not against the configured groups.* A k-mer
present in a retained-intron transcript of the same gene, or in a paralogue, is unique to
nothing as far as the quantifier is concerned. The gene's other transcripts are the default
background; an arbitrary FASTA — ideally the one the index was built from — can be streamed
in addition. Computed after the gate was graded, 1 of the 48 contrasts of the simulation
estimable as abundance against the gene's own transcripts was not estimable against the
reference-chromosome transcripts of its index, and 3 were not once that index's records on
scaffolds were added. The one against the reference chromosomes is SMN2's, whose summed TPM
did not pass the gate.

*The read model counts reads, not fragments.* A unique k-mer is informative only when it
falls inside a sequenced end; the expected count is evaluated exactly over start positions
and a truncated-normal fragment-length distribution.

*The compatibility system is a sequence-derived surrogate, not the observation model of a
sequencing run.* It collapses windows carried by the same set of transcripts at a window
length the user sets. We make no claim relating the two — not that a verdict here transfers
to reads, and not that it is conservative in either direction. The counterexample that ended
the last attempt to state such a relation is a regression test: four short transcripts whose
class contrast is estimable at window 3 and 4, *not* at 5 and 6, and again at 7 and 8. A
longer `window` is a different system, not a sharper one. Whether the surrogate predicts
what a quantifier recovers is an empirical question this paper does not answer. The reported
figure is the smallest |log2 fold change| a 95% interval excludes zero for, from the
delta-method SE of the log class ratio under the Poisson-weighted GLS covariance; the effect
size the user states, not the verdict, sets the exit status. A conditioning threshold,
`--tau`, was withdrawn rather than calibrated. Past `LINEARISATION_LIMIT`, a relative
standard error, a figure is flagged and counts as not resolved at any effect size.

Cohorts are combined by a weighted Stouffer combination [@stouffer1949; @liptak1958] of
per-cohort exact tests, the default, and by stratified signed-rank statistics with van
Elteren's weights [@vanelteren1960], beside the donor-pooled test, which lets depth or
tissue handling drive the result. Every test carries the finest p-value its design could
resolve, and every effect size a donor bootstrap interval.

A run can be pinned to an annotation release, its sequence saved and replayed offline.
Dependencies are the Ensembl REST API [@yates2015], NumPy [@harris2020], SciPy
[@virtanen2020] and Matplotlib [@hunter2007].

# Research impact statement

`isoform-dominance` is the generalised, reusable implementation of an isoform-usage
analysis developed for a study of leptin-receptor isoforms in the human choroid plexus.
The evidence offered here is concrete rather than aspirational.

*Use in research.* Version 2.1.1 performs the LEPR quantification and aggregation in
@yang2026, cited there in the Code Availability statement and the reference list, and
archived at a version DOI so the cited state stays reachable independently of this
package's later development. The author is a co-author of that study; the package is
otherwise distributed through PyPI.
<!-- If @yang2026 has not been accepted by submission, delete this paragraph and the
     bib entry: the section must stand on the evidence below without it. -->

*Reproducibility.* The bundled, download-free self-test regenerates the per-cohort
result that study reports — short-isoform predominance in control human choroid plexus
in each of two independent public cohorts, 5/5 and 6/6 donors [@gse228458; @gse137619] —
on a clean machine in seconds, and additionally reports the pooled, Stouffer and
stratified combinations across them, which that study does not (it analyses the cohorts
separately). The study's standalone reproducibility repository runs its analysis through
this package, so the result is regenerated by the software itself rather than by
separate one-off code.

<!-- AI usage disclosure: removed 091526 at author's request.
     Restore before submission from ../AI-DISCLOSURE-PARAGRAPH.md .
     JOSS requires the tools/models used and their versions. -->

# Acknowledgements

Computation used the O2 High Performance Compute Cluster at Harvard Medical School. The
author thanks the leptin-receptor study team for the data that motivated the bundled
example.

# References
