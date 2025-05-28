#!/bin/bash

################################################################################
# Script: 06-runAltAFplotter.sh
# Title: Run AltAFplotter Pipeline (nonSVs or SVs)
# Author: Marta Sevilla Porras
# Date: 26/05/2025
#
# Description:
# Runs UPD detection using ROH + inheritance pattern evaluation per individual
# using bcftools and a custom summarization script (`upd_finder.py`).
#
# Each simulation folder includes:
#  - One VCF for each of proband, father, and mother
#  - The aligned file list (aligned_files.csv)
#
# Requirements:
# - bcftools ≥ 1.12
# - Python ≥ 3.9 with:
#     pandas, altair, cyvcf2, multiprocess, openpyxl
# - upd_finder.py (from: https://github.com/HUGLeipzig/BatchUPDetection)
#
# Usage:
# sbatch 06-runAltAFplotter.sh
################################################################################

#SBATCH --array=0-4

# --------------------------
# Configuration
# --------------------------

REGION_SIZES=(1 2 5 10 20)
SIZE=${REGION_SIZES[$SLURM_ARRAY_TASK_ID]}

# Base paths
VCF_DIR="results_simulation/updbatch/${SIZE}mb/"
OUT_ROH="results_event/updbatch/${SIZE}mb/"
OUT_ISEC="results_event/updbatch/${SIZE}mb/"
OUT_RES="results_event/updbatch/${SIZE}mb/"

# CSV with list of aligned trio VCFs
FAMILY_FILE="${VCF_DIR}/aligned_files.csv"
UPD_SCRIPT="scripts/upd_finder.py"

mkdir -p "$OUT_ROH" "$OUT_ISEC" "$OUT_RES"

# --------------------------
# Process each trio
# --------------------------

header=$(head -n 1 "$FAMILY_FILE")

tail -n +2 "$FAMILY_FILE" | while IFS= read -r line; do
  [ -z "$line" ] && continue

  base_name=$(echo "$line" | cut -d',' -f1 | sed 's/_proband_merge_sort.vcf.gz//')
  echo "Processing: $base_name [${SIZE}mb]"

  out_folder="${OUT_RES}/${base_name}/"
  mkdir -p "$out_folder"

  temp_csv="${out_folder}/tmp.csv"
  echo "$header" > "$temp_csv"
  echo "$line" >> "$temp_csv"

  proband_vcf="${VCF_DIR}${base_name}_proband_merge_sort.vcf.gz"
  mother_vcf="${VCF_DIR}${base_name}_mother_merge_sort.vcf.gz"
  father_vcf="${VCF_DIR}${base_name}_father_merge_sort.vcf.gz"

  if [[ ! -f "$proband_vcf" || ! -f "$mother_vcf" || ! -f "$father_vcf" ]]; then
    echo "Missing VCFs for $base_name — skipping."
    continue
  fi

  roh_out="${OUT_ROH}/${base_name}_proband_roh.txt"
  isec_out="${OUT_ISEC}/${base_name}_proband_isec.txt"

  echo "Running bcftools roh..."
  bcftools roh -I "$proband_vcf" -G30 --AF-dflt 0.4 | awk '$1 == "RG" {print $0}' > "$roh_out"

  echo "Running bcftools isec..."
  bcftools isec -n +1 "$proband_vcf" "$mother_vcf" "$father_vcf" > "$isec_out"

  echo "Running Python summarization..."
  python "$UPD_SCRIPT" "$OUT_ROH" "$OUT_ISEC" "$temp_csv" "$VCF_DIR" "$out_folder"

  rm "$temp_csv"
  echo "Finished: $base_name"
done