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
receptor with a signalling-competent long isoform and a truncated short isoform that
retains ligand binding but loses the intracellular effector domain. Answering it from
bulk RNA-seq means assigning annotated transcripts to each class, quantifying each class
per sample, and testing the contrast across donors — and, before any of that,
establishing that short reads can resolve the comparison at all.

`isoform-dominance` packages these steps behind one Python command-line interface. Given
a gene symbol it proposes isoform groups by clustering the annotation's protein-coding
transcripts [@frankish2021] on their 3' terminal-exon splice acceptor. It then answers
the measurability question in three layers, before a single read is quantified: how much
sequence is unique to each class once the rest of the gene — or an arbitrary background
FASTA, such as the one a Salmon index was built from — is accounted for; how many
fragments would be sequenced across that unique sequence at a stated depth, read length
and fragment-length distribution, and hence what counting-noise floor the class ratio
carries; and whether each class total, and the contrast between them, is an *estimable*
function of the fragment-class system at all. Quantification is then read from `Salmon`
`quant.sf` files [@patro2017] and summarised per donor; the contrast is tested with a
donor-level Wilcoxon signed-rank test [@wilcoxon1945], reported per cohort with a
bootstrap interval, the finest p-value the design could have resolved, and a stratified
combination across cohorts. A bundled, download-free self-test reproduces the reference
result on a clean machine in seconds.

# Statement of need

The audience is bench and computational biologists testing a hypothesis about one gene
from bulk RNA-seq, who do not need a genome-wide screen but do need the focused
comparison done correctly and reproducibly.

Doing it by hand has a recurring, often-unstated pitfall. A short-read quantifier
apportions fragments among transcripts by solving a linear inverse problem, and a
transcript's abundance — or a group's — is recoverable only when it is an estimable
function of that system. An EM quantifier returns abundances for a non-identifiable
model as readily as for an identifiable one: the likelihood is flat along the
unidentified directions and the optimiser stops somewhere on that ridge. So a
fold-change reported between two classes the data cannot separate is not a weak result,
it is an artefact of tie-breaking — and it is indistinguishable in the output from a
real one. "The quantifier will work it out" is not an argument.

The obvious proxy — "does each class own at least one unique k-mer?" — is wrong in both
directions, and releases up to v2.1.1 used it. It is falsely reassuring: a class whose
only unique sequence is the thirty k-mers spanning one splice junction passes, while
yielding a handful of informative fragments at realistic depth. It is also falsely
alarming: a class nested inside another owns no unique k-mer at all, yet its total is
estimable, because the containing class is pinned by its own unique sequence and the
nested class's indicator therefore lies in the row space. Both failure modes are
reproduced as tests here.

`isoform-dominance` replaces the proxy with the estimability condition itself, evaluated
on the class-collapsed system and reported beside the two quantities that say whether an
estimable comparison is also a usable one: the conditioning of the contrast and the
expected number of informative fragments. It withholds a pass when the contrast lies
outside the row space, and says "weakly identifiable, and here is why" in the large
middle ground where the honest answer is neither yes nor no.

# State of the field

**Differential transcript usage.** DTU is a mature area whose tools target a different
task. `DEXSeq` [@anders2012] tests differential *exon* usage from exon-bin counts;
`DRIMSeq` [@nowicka2016], `satuRn` [@gilis2022] and `IsoformSwitchAnalyzeR`
[@vittingseerup2019] test transcript-level DTU genome wide; and `fishpond`/`swish`
[@zhu2019] propagates Salmon's quantification uncertainty into that testing. All of
them begin from a prepared count matrix, typically imported with `tximport`
[@soneson2016], with the transcript groups already defined by the analyst.

**Identifiability.** @hiller2009 gave criteria under which isoform abundances are
uniquely determined by junction-array and RNA-seq observations, and reported up to 97%
of 2,256 alternatively spliced human RefSeq genes identifiable from RNA-seq;
@ferrerbonsoms2022 extended these to a read- and fragment-length-aware criterion and
applied it genome wide. @zheng2022 measured what happens where they fail: for 20-47% of
Human Body Map transcripts the quantification error from nonidentifiability is large
enough that "the ranking of expression between the transcript and other isoforms from
the same gene cannot be determined" — precisely the comparison this package is asked to
certify. `terminus` [@sarkar2020] approaches the difficulty from the other side,
grouping transcripts *post hoc and from the data* by inferential uncertainty so that
analysis proceeds at whatever resolution the data support.

**The gap.** Those address *transcript-level* identifiability of the full
deconvolution, or the resolution the data happen to permit. Here the classes are fixed
*a priori and on biological grounds* — a signalling-competent receptor against a
truncated one — and the estimand is a contrast of class totals. Two questions then hide
under one word. A contrast is *estimable* when it lies in the row space of the
compatibility system: exact arithmetic, yes or no. It is *measurable* when it is
recoverable at finite sequencing depth. Transcript-level identifiability is not
necessary for the first — two transcripts of identical sequence are individually
unidentifiable while their sum is perfectly estimable, so a transcript-level verdict can
reject a comparison that is sound. Nor does it discriminate when it holds: full column
rank makes *every* contrast estimable, so the test stops separating cases exactly where
separation is needed. And it is silent on the second. That is why this package reports
an estimability verdict and a conditioning factor together, on the user's own contrast,
before quantification — and why neither is reported alone.

# Software design

The package is small, independently testable modules behind one `argparse` CLI:
`annotate` (Ensembl query and terminal-exon clustering), `identifiability` (the three
layers above), `extract` (Salmon `quant.sf` to per-donor group TPM), `stats` (signed-rank
test, cohort combination, plotting) and `contamination` (marker-correlation control).
Configuration is one JSON file naming the groups, so the pipeline applies to any gene
without code changes, and every subcommand accepts `--json`.

*Uniqueness is judged against a background, not against the configured groups.* A k-mer
absent from the other class but present in an unlisted retained-intron transcript of the
same gene, or in a paralogue, is unique to nothing as far as the quantifier is concerned.
The gene's other transcripts are fetched as background whenever the tool retrieves
sequence itself; supplied offline, background must be supplied with them. An arbitrary
FASTA — ideally the one the index was built from — is streamed without materialising its
own k-mer set. k-mers are folded to canonical form, because that is what the index
stores and what an unstranded library requires.

*The read model counts reads, not fragments.* A unique k-mer is informative only when it
falls entirely inside a sequenced end, so a stretch straddling the gap between a pair's
two reads, or lying beyond them in a longer fragment, contributes nothing. The expected
count is evaluated exactly over all start positions, averaged over a truncated-normal
fragment-length distribution, and converted to a counting-noise floor on the log2 class
ratio.

*The compatibility system is a sequence-derived surrogate, not the observation model of
a sequencing run.* It collapses windows carried by the same set of transcripts at a
window length the user sets; what a paired-end run observes is decided by read length,
the fragment-length distribution, and the fact that only the two ends of a fragment are
sequenced. We make no claim relating the two — not that a verdict here transfers to
reads, and not that it is conservative in either direction. Several attempts to state
such a relation failed in development; the counterexample that ended the last is now a
regression test, four short transcripts whose class contrast is estimable at window 3
and 4, *not* at 5 and 6, and estimable again at 7 and 8. A longer `window` is a different
system, not a sharper one. Whether the surrogate predicts what a quantifier recovers is
an empirical question this paper does not answer; calibration by simulation is tracked in
the repository issues, and the three verdict thresholds are provisional until it is done.
The conditioning factor `sqrt(c' (A'A)^+ c)` is reported with it: a standard-deviation
factor, not a variance, and the value it would take under `Var(y) = sigma^2 I`, which a
short-read quantifier does not satisfy. What it captures is geometry — a contrast
direction nearly degenerate in the observable classes is poorly determined however the
noise is distributed.

Cohorts are combined two ways beside the donor-pooled test earlier releases reported
alone: a weighted Stouffer combination [@stouffer1949; @liptak1958] of the per-cohort
exact tests, and a weighted combination of within-stratum signed-rank statistics using
van Elteren's `1/(n+1)` weights [@vanelteren1960] — the weighting only, since his test is
a stratified *two-sample* procedure and these data are paired within donor. Stouffer is
the default: at these sizes each stratum's exact p-value is trustworthy and the normal
approximation behind a combined rank statistic is not. Pooling donors across independent
studies ranks one cohort's differences against another's and lets depth or tissue
handling drive the result. Every per-cohort and pooled test carries the finest p-value
its design could resolve — under the exhaustive sign-permutation null one assignment puts
every difference on the same side, so a two-sided test on n non-zero pairs cannot report
below 2^(1-n), which at n = 5 is 0.0625 — and every effect size a donor bootstrap
interval.

`identifiability` accepts transcript sequences and background offline, keeping the test
suite deterministic and allowing use on compute nodes without internet access; `annotate`
queries Ensembl and has no offline mode. Dependencies are minimal: the Ensembl REST API
[@yates2015], NumPy [@harris2020], SciPy [@virtanen2020] and Matplotlib [@hunter2007].

# Research impact

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

*Generality.* The method is demonstrated on three biologically unrelated loci — the
*NTRK2*, *NTRK3* and *FLT1* worked examples in the documentation — each confirming the
proposed classes are short-read separable.

# AI usage disclosure

Generative AI assistance (Anthropic Claude; **[CONFIRM BEFORE SUBMISSION: exact model
names and versions used for v1.0–v2.1]**, and Claude Opus 5 for the v2.2 identifiability
and statistics work) was used for code scaffolding, refactoring, test drafting,
documentation, and copy-editing of this paper. All AI-assisted outputs were reviewed,
edited, and validated by the author, who made the core design decisions — the
terminal-exon grouping approach; the decision to evaluate estimability of the
class-collapsed system rather than transcript-level identifiability; the choice of
conditioning and expected-informative-read thresholds, and their provisional status
pending simulation; the read-level rather than
fragment-level informativeness model; the stratified cohort combination; and the
package architecture — and is fully responsible for the correctness, originality, and
licensing of the submitted materials.

# Acknowledgements

Computation used the O2 High Performance Compute Cluster at Harvard Medical School. The
author thanks the leptin-receptor study team for the data that motivated the bundled
example.

# References
