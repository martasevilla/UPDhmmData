#!/bin/bash

################################################################################
# Script: 05-runUPDhmm.sh
# Title: Run UPDhmm on Simulated VCFs (nonSVs or SVs)
# Author: Marta Sevilla Porras
# Date: 26/05/2025
#
# Description:
# This script runs the `calculateEvents()` function from the UPDhmm R package
# on VCFs simulated in step 03.
#
# Requirements:
# - R ≥ 4.2 with packages: UPDhmm, VariantAnnotation
#
# Inputs:
# - Simulated VCFs: results_simulation/updhmm/<size>/
#
# Outputs:
# - CSV files with UPD event predictions: results_event/updhmm/<size>/
#
# Usage:
# bash 05-runUPDhmm.sh
################################################################################

# --------------------------
# Configuration
# --------------------------


INPUT_BASE="results_simulation/updhmm"
OUTPUT_BASE="results_event/updhmm"
REGION_SIZES=("1mb" "2mb" "5mb" "10mb" "20mb")

module load R/4.2.0-foss-2021b

# --------------------------
# Process each region size
# --------------------------
for size in "${REGION_SIZES[@]}"; do

  input_dir="${INPUT_BASE}/${size}"
  output_dir="${OUTPUT_BASE}/${size}"
  mkdir -p "$output_dir"

  echo "Processing size: $size"

  shopt -s nullglob
  for vcf_file in "$input_dir"/*merge_sort.vcf.gz; do
    file_name=$(basename "$vcf_file")
    csv_file="${output_dir}/events_${file_name%.vcf.gz}.csv"

    echo "Running UPDhmm on: $file_name"

    Rscript -e "
      library(UPDhmm)
      library(VariantAnnotation)
      vcf <- readVcf('$vcf_file')
      colnames(vcf) <- c('proband', 'father', 'mother')
      events <- calculateEvents(vcf)
      write.csv(events, file = '$csv_file', row.names = FALSE)
    "

    echo "Output saved to: $csv_file"
  done

done