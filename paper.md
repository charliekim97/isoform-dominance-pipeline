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

Many biological questions reduce to a single, focused comparison: *for one gene, which
functional isoform class predominates in a given tissue or condition?* A common example
is a receptor with a signalling-competent long isoform and a truncated short isoform
that retains ligand binding but loses the intracellular effector domain. Answering this
from bulk RNA-seq requires deciding which annotated transcripts belong to each
functional class, quantifying each class per sample, and testing the contrast across
donors — and, before any of that, establishing that short reads can resolve the
comparison at all.

`isoform-dominance` packages these steps behind a single Python command-line interface.
Given a gene symbol it proposes isoform groups from Ensembl by clustering
protein-coding transcripts on their 3' terminal-exon splice acceptor [@frankish2021].
It then answers the measurability question in three layers, all computed from sequence
alone and before a single read is quantified: how much sequence is unique to each class
once the rest of the gene — or an arbitrary background FASTA, such as the one a Salmon
index was built from — is taken into account; how many fragments would actually be
sequenced across that unique sequence at a stated depth, read length and
fragment-length distribution, and hence what counting-noise floor the class ratio
carries; and whether each class total, and the contrast between the two class totals,
is an *estimable* function of the fragment-class system at all. Quantification is then
read from `Salmon` `quant.sf` files [@patro2017] and summarised per donor; the class
contrast is tested with a donor-level exact Wilcoxon signed-rank test [@virtanen2020],
reported per cohort with a bootstrap interval on the effect size, the finest p-value the
design could have resolved, and a stratified combination across cohorts. A bundled,
download-free self-test reproduces the reference result on a clean machine in seconds.

# Statement of need

The target audience is bench and computational biologists testing a hypothesis about
one gene from bulk RNA-seq, who do not need a genome-wide screen but do need the
focused comparison done correctly and reproducibly.

Doing this by hand has a recurring, often-unstated pitfall. A short-read quantifier
apportions fragments among transcripts by solving a linear inverse problem, and the
abundance of a transcript — or of a group of transcripts — is recoverable only when it
is an estimable function of that system. A fold-change reported between two isoform
classes that the data cannot separate is not a weak result; it is an artefact of the
optimiser's tie-breaking, and it is indistinguishable in the output from a real one.

The obvious proxy for this — "does each class own at least one unique k-mer?" — is
wrong in both directions, and releases of this package up to v2.1.1 used it. It is
falsely reassuring: a class whose only unique sequence is the thirty k-mers spanning
one splice junction passes, while yielding a handful of informative fragments at
realistic depth. It is also falsely alarming: a class nested inside another owns no
unique k-mer at all, yet its total is still an estimable function of the fragment-class
system, because the containing class is pinned by its own unique sequence and the
indicator of the nested class therefore lies in the row space. Both failure modes are
reproduced as tests in this package.

The distinction matters for how the verdict is justified. An EM quantifier returns
abundances for a non-identifiable model as readily as for an identifiable one -- the
likelihood is flat along the unidentified directions and the optimiser stops somewhere
on that ridge -- so "the quantifier will work it out" is not an argument that a class
total is recoverable. What settles it is whether the coefficient vector lies in the row
space of the design matrix, and how well conditioned that direction is; that is
what this package computes and reports.

`isoform-dominance` replaces that proxy with the estimability condition itself,
evaluated on the class-collapsed system, and reports it beside the two quantities that
indicate whether an estimable comparison is also a usable one: the conditioning of the
contrast and the expected number of informative fragments. It withholds a pass when the
contrast lies outside the row space of that system, and says "weakly identifiable, and
here is why" in the large middle ground where the honest answer is neither yes nor no.

# State of the field

**Differential transcript usage.** DTU is a mature area with strong tools that target a
different task and audience. `DEXSeq` [@anders2012], `DRIMSeq` [@nowicka2016] and
`satuRn` [@gilis2022] perform genome-wide DTU testing and assume the user has already
produced a transcript-by-sample count matrix and defined the transcript groups.
`IsoformSwitchAnalyzeR` [@vittingseerup2019] adds rich functional annotation of isoform
switches but is an R/Bioconductor workflow in which transcript grouping and import are
configured by the analyst. `fishpond`/`swish` [@zhu2019] rigorously propagates
quantification uncertainty using Salmon inferential replicates, but likewise operates
on a prepared dataset. Transcript-level import for these workflows is typically handled
by `tximport` [@soneson2016].

**Identifiability.** The identifiability of the isoform deconvolution problem is itself
a studied question, and this package builds on that literature rather than
rediscovering it. @hiller2009 gave criteria under which isoform abundances are uniquely
determined by junction-array and RNA-seq observations, and reported that most
alternatively spliced human genes are identifiable from RNA-seq. @ferrerbonsoms2022
extended this to a read-length- and fragment-length-aware criterion, applied it genome
wide, and used it to choose library fragment lengths. Both address *transcript-level*
identifiability of the full deconvolution. `terminus` [@sarkar2020] approaches the same
difficulty from the other side, grouping transcripts *post hoc and from the data* by
inferential uncertainty, so that analysis proceeds at whatever resolution the data
support.

**Build vs. contribute.** The gap those leave is the one a biologist actually stands
in. The classes here are defined *a priori and on biological grounds* — a
signalling-competent receptor against a truncated one — and the estimand is the class
total and the contrast between two class totals, not the individual transcript
abundances. Two questions hide under one word here, and separating them is the point. A
contrast is *estimable* when it lies in the row space of the compatibility system — a
yes/no in exact arithmetic. It is *measurable* when it can be recovered at a finite
sequencing depth. Transcript-level identifiability settles the first in one direction
only, and is silent on the second.

It is not necessary for estimability: two transcripts with identical sequence make each
abundance unidentifiable while their sum is perfectly estimable, so a transcript-level
verdict can reject a comparison that is in fact sound.

It is sufficient for estimability, and that is the problem rather than the reassurance
it sounds like. If every transcript is identifiable the design matrix has full column
rank, so *every* contrast is estimable; the row-space condition is satisfied trivially
and a rank verdict stops discriminating exactly where it is needed. It is not sufficient
for measurability: full rank carries no statement about finite depth, and the class
direction can be so weakly observed that a formally estimable contrast is unrecoverable
in practice. That is why the estimability verdict is reported together with a
conditioning factor, and why neither is reported alone.

`terminus` will not answer the question because it chooses the groups itself, after
seeing the data, and so cannot say whether *the comparison the investigator came with*
is supported; the transcript-level criteria will not answer it because their failing
verdict does not transfer — a gene they declare non-identifiable may still determine the
class contrast exactly. `isoform-dominance` evaluates the estimability of the user's own
class contrast, before quantification, and none of the DTU tools above takes a bare
gene symbol as input or performs such a check. It is deliberately narrow and
complementary rather than competing: it does not attempt genome-wide discovery,
delegates quantification to Salmon rather than re-implementing it, and presents group
proposals for review rather than as final calls.

The tool generalises across unrelated biology. Beyond the bundled leptin-receptor
example, `annotate` + `identifiability` cleanly recover, and confirm as short-read
separable: the full-length vs truncated kinase receptors of *NTRK2* (TrkB-FL vs
TrkB.T1) and *NTRK3* (TrkC), and the membrane vs secreted-decoy forms of *FLT1* (full
VEGFR1 vs soluble sFlt-1) — worked examples included in the documentation.

# Software design

The package is organised as small, independently testable modules behind a single
`argparse` CLI: `annotate` (Ensembl REST query and terminal-exon clustering),
`identifiability` (the three layers described above), `extract` (Salmon `quant.sf`
parsing and per-donor group TPM), `stats` (paired exact Wilcoxon, cohort combination
and plotting), and `contamination` (marker-correlation control). Configuration is a
single JSON file describing the isoform groups, so the same pipeline applies to any
gene without code changes, and every subcommand accepts `--json` so the tool composes
inside a workflow manager without anyone parsing its printed output.

Three design choices are worth stating because they are where the package's judgement
sits rather than its plumbing.

*Uniqueness is judged against a background, not against the configured groups.* A
k-mer absent from the other class but present in an unlisted retained-intron transcript
of the same gene, or in a paralogue, is not unique to anything as far as the quantifier
is concerned. The gene's remaining transcripts are used by default, and an arbitrary
FASTA — ideally the one the index was built from — can be streamed as background in
memory proportional to the query rather than the file. k-mers are folded to their
canonical form, because that is what the index stores and what an unstranded library
requires.

*The read model counts reads, not fragments.* A unique stretch is informative only if a
sequenced end covers it; with 200 nt fragments and 100 nt reads a unique region in the
middle of the fragment is never observed. The expected informative-fragment count is
evaluated exactly over all start positions and averaged over a truncated-normal
fragment-length distribution, then converted to a counting-noise floor on the log2
class ratio.

*The compatibility system is a sequence-derived surrogate, not the observation model
of a sequencing run.* It collapses windows carried by the same set of transcripts, at a
window length the user sets; what an actual paired-end run observes is decided by read
length, the fragment-length distribution, and the fact that only the two ends of a
fragment are sequenced.

We make no claim relating the two — not that a verdict here transfers to reads, and not
that it is conservative in either direction. Several attempts to state such a relation
failed in development, and the counterexample that ended the last of them is now a
regression test: four short transcripts whose class contrast is estimable at window 3
and 4, *not* estimable at 5 and 6, and estimable again at 7 and 8. A longer `window` is
a different system, not a sharper one.

What the verdict is, exactly: `c` lies in the row space of a compatibility system built
from sequence at a stated window length, reported with a structural conditioning factor
because a full-rank system with a near-degenerate contrast direction passes a rank test
and still yields nothing. Whether that predicts what a quantifier recovers is an
empirical question this paper does not answer; calibration by simulation is planned
work, tracked in the repository issues.

That conditioning factor, `sqrt(c' (A'A)^+ c)`, is deliberately not called a variance —
note the square root: the variance factor is `c'(A'A)^+ c` and this is its
standard-deviation counterpart. It is the value that standard-deviation factor would
take under `Var(y) = sigma^2 I`, and a short-read quantifier does not satisfy that: fragment counts
are heteroskedastic, and Salmon's rich equivalence classes carry per-transcript
conditional probabilities and bias weights rather than the compatibility-only
construction used here. What the quantity captures is geometry — a contrast direction
that is nearly degenerate in the observable classes is poorly determined however the
noise is distributed. Whether it predicts realised quantification error is an empirical
question, addressed by simulation rather than asserted here, and the thresholds
separating the three verdicts are provisional until that calibration is done.

On the statistics, cohorts are combined two ways beside the donor-pooled test that
earlier releases reported alone: a weighted Stouffer combination of the per-cohort exact
tests [@stouffer1949], and a weighted combination of the within-stratum signed-rank
statistics using van Elteren's design-free `1/(n+1)` weights [@vanelteren1960] -- the
weighting scheme only, since van Elteren's test is a stratified *two-sample* procedure
and these data are paired within donor. Stouffer is the default of the two, because at
these sample sizes each stratum's exact p-value is trustworthy and the normal
approximation behind a combined rank statistic is not. Pooling donors from independent
studies ranks one cohort's differences against another's and lets depth or tissue
handling drive the result. Every test is reported with the finest p-value its own design could have resolved: under
the exhaustive sign-permutation null exactly one assignment puts every difference on the
same side, so a two-sided test on n non-zero pairs cannot report below 2^(1-n), and at
n = 5 that is 0.0625. Ties among the absolute differences do not raise this bound. Every
effect size is reported with a donor bootstrap interval.

Network-dependent steps accept offline inputs to keep the test suite deterministic and
to allow use on compute nodes without internet access. Dependencies are kept minimal:
the annotation and sequence queries use the Ensembl REST API [@yates2015], numerical
work uses NumPy [@harris2020] and SciPy [@virtanen2020], and figures use Matplotlib
[@hunter2007].

# Research impact

`isoform-dominance` is the generalised, reusable implementation of an isoform-usage
analysis developed for a study of leptin-receptor isoforms in the human choroid plexus,
and the evidence offered here is concrete rather than aspirational.

First, use in research: version 2.1.1 performs the LEPR isoform quantification and
aggregation in @yang2026, cited there both in the Code Availability statement and in the
reference list, and archived at a version DOI so that the cited state remains reachable
independently of this package's later development. The author of this package is a
co-author of that study; the package is otherwise distributed through PyPI.
<!-- If @yang2026 has not been accepted by submission, delete this paragraph and the
     bib entry: the section must stand on the reproducibility and benchmark evidence
     below without it. -->

Second, reproducibility: the bundled, download-free self-test regenerates the
per-cohort result that study reports — short-isoform predominance in control human
choroid plexus in each of two independent public cohorts, 5/5 and 6/6 donors
[@gse228458; @gse137619] — on a clean machine in seconds, and additionally reports the
pooled and stratified combinations across those cohorts, which that study does not
(it analyses the cohorts separately). The study's standalone reproducibility repository
runs its analysis through this package, so the result is regenerated by the software
itself rather than by separate one-off code.

Third, generality: the method is demonstrated on three further, biologically unrelated
loci with verifiable outputs — the *NTRK2*, *NTRK3* and *FLT1* worked examples in the
documentation — each confirming the proposed isoform classes are short-read separable.

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
