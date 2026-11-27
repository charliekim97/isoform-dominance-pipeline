# Assessing isoform-class measurements from short-read RNA-seq

Short reads can leave individual transcripts indistinguishable while still identifying the
total abundance of a biologically defined class. `isoform-dominance` uses a sequence-derived
observation model to assess which class totals and contrasts are structurally estimable and
to predict counting precision for a specified read design and depth. The workflow also
aggregates Salmon transcript estimates and tests class differences across donors.

## What the benchmarks showed

**Coverage skew can bias a structurally estimable contrast.** The simulation used 49 genes
and seven positional coverage designs. Salmon and kallisto gave evaluable results in every
design for 47 genes. Among the 39 estimable Salmon fragment-count contrasts with finite
predicted standard errors, absolute bias exceeded the predicted standard error in **0 of 39**
under uniform coverage and **20 to 24 of 39** under mild positional skew, depending on
its direction. Replicate standard deviation stayed below the prediction in 34 to 37 of 39
under the same mild skew. [Simulation results](validation/README.md#results)

**SMN2 stopped a proposed default change.** A prespecified test of the length-corrected
abundance estimand used class ratios from summed TPM. Six of seven primary contrasts met
both limits for bias and replicate variation. SMN2 exceeded both, so the fragment estimand
stayed the default. A follow-up held the simulated libraries fixed and removed 33 SMN1 and
SMN2 alternate-locus records from the index while retaining the chromosome 5 SMN1 records.
SMN2's absolute log2-ratio bias fell from **0.474 to 0.110** and its replicate standard
deviation from **0.751 to 0.352**, both below the prespecified standard error of 0.442.
These libraries contained no SMN1-derived fragments. [Gate and follow-up](validation/README.md#d1s-result)

## Use in research

Version 2.1.1 was used for the LEPR analysis in Yang et al. The offline self-test replays
aggregation and statistics from stored transcript TPM, reproducing short-class predominance
in 5 of 5 and 6 of 6 control donors in two human choroid-plexus cohorts.
[Study](https://doi.org/10.1038/s41467-026-79556-9) ·
[Archived software](https://doi.org/10.5281/zenodo.20738150) ·
[Raw-read workflow](https://github.com/charliekim97/Human-choroid-plexus-LEPR-isoform-analysis-reproducible-code)

The benchmark covers a fixed depth and expression scheme in simulated libraries.
Measuring accuracy in biological samples requires independent quantitative evidence.
