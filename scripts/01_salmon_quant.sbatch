#!/bin/bash
#SBATCH --job-name=isoform_quant
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=isoform_quant_%j.out
#SBATCH --error=isoform_quant_%j.err
# =====================================================================
# 01 - Salmon quantification (generic; SLURM/HPC, e.g. Harvard O2).
# Downloads paired FASTQ from ENA and runs salmon quant per donor.
# Reuses a pre-built decoy-aware index (build once, see README).
# =====================================================================
set -euo pipefail

# ============ EDIT THESE ============
SALMON=salmon                              # full path to salmon binary if not on PATH
INDEX=/n/scratch/.../decoy_index           # pre-built decoy-aware index (see README)
SAMPLE_MAP=sample_map.csv                   # columns: donor,condition,SRR
OUTDIR=$PWD/quant
# ===================================
THREADS=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$OUTDIR" fastq
command -v "$SALMON" >/dev/null || { echo "ERROR: salmon not found ($SALMON)"; exit 1; }
"$SALMON" --version

ena(){ local s=$1; echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${s:0:6}/0${s: -2}/${s}"; }
mapfile -t LINES < <(tail -n +2 "$SAMPLE_MAP")

# --- Phase 1: parallel FASTQ download (wget -c resumes/skips complete) ---
: > urls.txt
for ln in "${LINES[@]}"; do
  IFS=, read -r donor cond srr <<< "$ln"
  b=$(ena "$srr"); echo "$b/${srr}_1.fastq.gz"; echo "$b/${srr}_2.fastq.gz"
done > urls.txt
echo ">>> downloading $(wc -l < urls.txt) files (6 parallel)"
xargs -P 6 -n 1 wget -c -q -P fastq/ < urls.txt || true

# --- Phase 2: quant each donor (idempotent; single-end R2 fallback) ---
for ln in "${LINES[@]}"; do
  IFS=, read -r donor cond srr <<< "$ln"
  if [ -f "$OUTDIR/$donor/quant.sf" ]; then echo ">>> $donor done, skip"; continue; fi
  echo ">>> quant $donor ($srr)"
  if [ -s fastq/${srr}_1.fastq.gz ] && [ -s fastq/${srr}_2.fastq.gz ]; then
    "$SALMON" quant -i "$INDEX" -l A -1 fastq/${srr}_1.fastq.gz -2 fastq/${srr}_2.fastq.gz \
      -p "$THREADS" --validateMappings --gcBias -o "$OUTDIR/$donor"
  elif [ -s fastq/${srr}_2.fastq.gz ]; then
    "$SALMON" quant -i "$INDEX" -l A -r fastq/${srr}_2.fastq.gz \
      -p "$THREADS" --validateMappings -o "$OUTDIR/$donor"
  else
    echo "!!! $donor ($srr): no fastq, skipped"
  fi
done
echo "DONE -> $OUTDIR/<donor>/quant.sf"
