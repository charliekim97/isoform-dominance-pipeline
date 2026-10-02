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
# ===== EDIT (or pass with sbatch --export=ALL,SAMPLE_MAP=...,OUTDIR=...) =====
SALMON=${SALMON:-salmon}
INDEX=${INDEX:-/path/to/decoy_index}
SAMPLE_MAP=${SAMPLE_MAP:-sample_map.csv}          # donor,condition,SRR
# one directory per sample map: two cohorts may reuse donor names (both examples have
# ctrl1-ctrl5), and a donor whose quant.sf exists is skipped below
OUTDIR=${OUTDIR:-$PWD/quant/$(basename "$SAMPLE_MAP" .csv)}
# ================
THREADS=${SLURM_CPUS_PER_TASK:-12}; mkdir -p "$OUTDIR" fastq
# ENA's directory for a run: vol1/fastq/<first 6>/<by run-number length>/<run>. The second
# level is absent for a 6-digit number, 00+last digit for 7, 0+last two for 8, the last
# three for 9 (ERR2060213 -> ERR206/003). bash 3.2 syntax: tests/test_sbatch.py runs it.
ena() {
  local s=$1 n=${#1} d
  n=$((n - 3))
  case $n in
    6) d="" ;;
    7) d="00${s: -1}/" ;;
    8) d="0${s: -2}/" ;;
    9) d="${s: -3}/" ;;
    *) echo "ena: cannot place run $s (expected 3 letters and 6-9 digits)" >&2; return 1 ;;
  esac
  echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${s:0:6}/${d}${s}"
}
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
