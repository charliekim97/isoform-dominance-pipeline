#!/bin/bash
#SBATCH --job-name=isoform_quant
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=isoform_quant_%j.out
# Salmon quantification (SLURM/HPC). Reuses a pre-built decoy-aware index.
# Build that index
#  - from the GENCODE release matching `annotate --ensembl-release`
#    (https://www.gencodegenes.org/human/releases.html), and
#  - from reference-chromosome transcripts only. GENCODE >= 48 transcripts.fa.gz and
#    Ensembl cdna.all.fa.gz (release 116 checked) also hold scaffold, patch and
#    alternate-locus copies, which take reads that `extract` does not count; filter
#    them out first (README, "Index scope"). Every donor of a cohort must be quantified
#    against the same index: `extract` compares aux_info/meta_info.json across donors
#    and stops if they differ.
set -euo pipefail
# ===== EDIT =====
SALMON=salmon
INDEX=/path/to/decoy_index
SAMPLE_MAP=sample_map.csv          # donor,condition,SRR
OUTDIR=$PWD/quant
# ================
THREADS=${SLURM_CPUS_PER_TASK:-12}; mkdir -p "$OUTDIR" fastq
ena(){ local s=$1; echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${s:0:6}/0${s: -2}/${s}"; }
mapfile -t LINES < <(tail -n +2 "$SAMPLE_MAP")
: > urls.txt
for ln in "${LINES[@]}"; do IFS=, read -r d c srr <<< "$ln"; b=$(ena "$srr"); echo "$b/${srr}_1.fastq.gz"; echo "$b/${srr}_2.fastq.gz"; done > urls.txt
# Download all FASTQs; fail the job if any download fails (reproducibility > silent skip).
xargs -P 6 -n 1 wget -c -q -P fastq/ < urls.txt
for ln in "${LINES[@]}"; do
  IFS=, read -r d c srr <<< "$ln"
  [ -f "$OUTDIR/$d/quant.sf" ] && continue
  for mate in 1 2; do
    [ -s "fastq/${srr}_${mate}.fastq.gz" ] || { echo "ERROR: missing or empty fastq/${srr}_${mate}.fastq.gz" >&2; exit 1; }
  done
  "$SALMON" quant -i "$INDEX" -l A -1 fastq/${srr}_1.fastq.gz -2 fastq/${srr}_2.fastq.gz \
    -p "$THREADS" --validateMappings --gcBias -o "$OUTDIR/$d"
done
echo "DONE -> $OUTDIR/<donor>/quant.sf"
