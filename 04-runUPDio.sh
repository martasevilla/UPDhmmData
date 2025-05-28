#!/bin/bash

################################################################################
# Script: 04-runUPDio.sh
# Title: Run UPDio on Simulated VCFs (nonSVs or SVs)
# Author: Marta Sevilla Porras
# Date: 26/05/2025
#
# Description:
# This script runs the UPDio algorithm on simulated multi-sample VCFs
# produced in step 03. It processes all region sizes (1–20Mb) for
# VCFs with SV filtering ("SVs")
#
# The UPDio Perl script is called using sample names extracted
# directly from the VCF header using bcftools.
#
# Requirements:
# - bcftools ≥ 1.16
# - Perl script UPDio.pl (https://github.com/findingdan/UPDio)
#
# Inputs:
# - Simulated VCFs: results_simulation/updio/<size>/
# - UPDio script:   path to UPDio.pl
#
# Output:
# - UPDio result files saved alongside each input VCF
#
# Usage:
# bash 04-runUPDio.sh
################################################################################

# --------------------------
# Configuration
# --------------------------

# Set mode: "nonSVs" or "SVs"


BASE_DIR="results_simulation/updio"
UPDIO_SCRIPT="/path/to/UPDio.pl"  # ← Replace with actual path to the script

# Region sizes to process
SIZES=("1mb" "2mb" "5mb" "10mb" "20mb")

# --------------------------
# Run UPDio on each VCF
# --------------------------
for size in "${SIZES[@]}"; do

  echo "Processing size: $size"

  INPUT_DIR="${BASE_DIR}/${size}"
  OUTPUT_DIR="${INPUT_DIR}"

  shopt -s nullglob
  for vcf_file in "$INPUT_DIR"/*merge_sort.vcf.gz; do

    echo "$(basename "$vcf_file")"

    # Extract trio sample names
    sample_names=($(bcftools query -l "$vcf_file"))
    if [[ ${#sample_names[@]} -lt 3 ]]; then
      echo "Error: Less than 3 samples in $vcf_file"
      continue
    fi

    proband="${sample_names[0]}"
    father="${sample_names[1]}"
    mother="${sample_names[2]}"

    # Run UPDio
    perl "$UPDIO_SCRIPT" \
      --multisample_vcf "$vcf_file" \
      --output_path "$OUTPUT_DIR" \
      --childID "$proband" \
      --momID "$mother" \
      --dadID "$father" \
      --name "$(basename "$vcf_file")" \
      --significance_level 0.05

    echo "Done: $vcf_file"

  done
  shopt -u nullglob

done