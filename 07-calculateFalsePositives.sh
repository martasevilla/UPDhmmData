#!/bin/bash

################################################################################
# Script: 08-calculateFalsePositives_cleanData.sh
# Title: Detect False Positives in Clean Trio Data (No Simulated UPDs)
# Author: Marta Sevilla Porras
# Date: 28/05/2025
#
# Description:
# Applies UPDhmm, UPDio and AltAFplotter to original, non-simulated VCFs
# to estimate the number of false positive UPD detections.
#
# Input:
#   - data/01_centro_seg_dup/*vcf.gz         (region-only filtered)
#   - data/01_centro_seg_dup_SV/*vcf.gz      (region + SV filtered)
#
# Output:
#   - Number of false positives detected by each method
#
# Usage:
#   bash 08-calculateFalsePositives_cleanData.sh data/01_centro_seg_dup/ or data/01_centro_seg_dup_SV/
################################################################################

INPUT_DIR="$1"
OUT_DIR="$1"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Directory '$INPUT_DIR' not found"
  echo "Usage: bash $0 <vcf_directory>"
  exit 1
fi

echo "=== FALSE POSITIVES FROM CLEAN VCFs ==="
echo "Input: $INPUT_DIR"
echo ""

mkdir -p "${OUT_DIR}/results_event/fp_updhmm" "${OUT_DIR}/results_event/fp_updio" "${OUT_DIR}/results_event/fp_altaf"

##########################################
# 1. UPDhmm
##########################################
echo "→ Running UPDhmm..."

for vcf_file in ${INPUT_DIR}/*.vcf.gz; do
  file_name=$(basename "$vcf_file")
  csv_file="${OUT_DIR}/results_event/fp_updhmm/${file_name%.vcf.gz}.csv"

  echo "Running UPDhmm on: $file_name"

  Rscript -e "
    library(UPDhmm)
    library(VariantAnnotation)
    vcf <- readVcf('$vcf_file')
    colnames(vcf) <- c('proband', 'father', 'mother')
    events <- calculateEvents(vcf)
    write.csv(events, file = '$csv_file', row.names = FALSE)
  "

done

##########################################
# 2. UPDio
##########################################
echo "→ Running UPDio..."

for vcf_file in ${INPUT_DIR}/*.vcf.gz; do
  file_name=$(basename "$vcf_file")
  samples=($(bcftools query -l "$vcf_file"))
  proband="${samples[0]}"
  father="${samples[1]}"
  mother="${samples[2]}"

  perl /path/to/UPDio.pl \
    --multisample_vcf "$vcf_file" \
    --output_path "${OUT_DIR}/results_event/fp_updio/" \
    --childID "$proband" --dadID "$father" --momID "$mother" \
    --name "${file_name%.vcf.gz}" --significance_level 0.05

done

##########################################
# 3. AltAFplotter
##########################################
echo "→ Running AltAFplotter..."

for vcf_file in ${INPUT_DIR}/*.vcf.gz; do
  file_name=$(basename "$vcf_file")
  base_name="${file_name%.vcf.gz}"
  samples=($(bcftools query -l "$vcf_file"))
  proband="${samples[0]}"
  father="${samples[1]}"
  mother="${samples[2]}"

  out_folder="${OUT_DIR}/results_event/fp_altaf/${base_name}/"
  mkdir -p "$out_folder"

  # Create individual VCFs
  for i in {0..2}; do
    role=("proband" "father" "mother")
    sample="${samples[$i]}"
    label="${role[$i]}"
    indiv_vcf="${out_folder}/${base_name}_${label}_prov.vcf.gz"
    final_vcf="${out_folder}/${base_name}_${label}_merge_sort.vcf.gz"

    bcftools view -s "$sample" "$vcf_file" -Oz -o "$indiv_vcf"
    tabix "$indiv_vcf"

    bcftools view -e 'FORMAT/GT=="0|0"' "$indiv_vcf" -Oz -o "$final_vcf"
    tabix "$final_vcf"
    rm "$indiv_vcf"*
  done

  roh_out="${out_folder}/${base_name}_proband_roh.txt"
  isec_out="${out_folder}/${base_name}_proband_isec.txt"

  echo "Running bcftools roh..."
  bcftools roh -I "${out_folder}/${base_name}_proband_merge_sort.vcf.gz" -G30 --AF-dflt 0.4 | awk '$1 == "RG" {print $0}' > "$roh_out"

  echo "Running bcftools isec..."
  bcftools isec -n +1 "${out_folder}/${base_name}_proband_merge_sort.vcf.gz" \
                     "${out_folder}/${base_name}_mother_merge_sort.vcf.gz" \
                     "${out_folder}/${base_name}_father_merge_sort.vcf.gz" > "$isec_out"

  echo "Running Python summarization..."
  python scripts/upd_finder.py "$out_folder" "$roh_out" "$isec_out" "$vcf_file" "$out_folder"

done