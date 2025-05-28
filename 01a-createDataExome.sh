#!/bin/bash

################################################################################
# Script: 01-createData.sh
# Title: Preprocessing Exome VCFs for Trio-Based UPD Simulations
# Author: Marta Sevilla Porras
# Date: 26/05/2025
#
# Description:
# This script preprocesses exome VCFs for selected trios by:
#  1. Downloading per-chromosome VCFs for each trio from UCSC
#  2. Creating complete multi-chromosome VCFs per trio
#  3. Generating, in parallel:
#     - VCF files ready for UPD simulation (EXOME, biallelic, region filtered)
#     - TSV files with structural variants (SVs) to be excluded
#
# Outputs:
#  - Filtered VCFs:
#    - data/01_centro_seg_dup/         (region-only filtered)
#    - data/01_centro_seg_dup_SV/      (region + SV filtered)
#  - TSV files: SVs_coordinates/<trio>_SV_filtered_not_INV.tsv
#
# Requirements:
#  - bcftools ≥ 1.19
#  - BED files: regions_to_exclude.bed
################################################################################

# -----------------------------
# Configuration
# -----------------------------
TRIOS=("NA19675" "NA19685")
CHROMOSOMES=$(seq 1 22)
UCSC_BASE_URL="https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes"
BED_REGIONS="regions_to_exclude.bed"
OUT_DIR_RAW="data/00_raw"
OUT_DIR_REGIONS="data/01_centro_seg_dup" #see supplementary material
OUT_DIR_SV="data/01_centro_seg_dup_SV"
SV_COORD_DIR="SVs_coordinates"

mkdir -p "$OUT_DIR_RAW" "$OUT_DIR_REGIONS" "$OUT_DIR_SV" "$SV_COORD_DIR"

# -----------------------------
# Process each trio
# -----------------------------
for trio in "${TRIOS[@]}"; do

  echo "Processing trio: $trio"

  # Define trio samples manually (update if needed)
  if [[ "$trio" == "NA19675" ]]; then
    sample_list="NA19675,NA19679,NA19678"  # proband,father,mother
    proband="NA19675"
  elif [[ "$trio" == "NA19685" ]]; then
    sample_list="NA19685,NA19661,NA19660"  # proband,father,mother
    proband="NA19685"
  fi

  # --------------------------------------------
  # Step 1: Download per-chromosome VCFs from UCSC and extract samples
  # --------------------------------------------
  for chr in $CHROMOSOMES; do
    bcftools view \
      -s $sample_list \
      -v snps,indels \
      --force-samples \
      -Oz -o "${OUT_DIR_RAW}/${trio}_chr${chr}.vcf.gz" \
      "${UCSC_BASE_URL}/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
  done

  # --------------------------------------------
  # Step 2: Concatenate per-chromosome VCFs into one file per trio
  # --------------------------------------------
  ls ${OUT_DIR_RAW}/${trio}_chr*.vcf.gz > ${OUT_DIR_RAW}/${trio}_file_list.txt
  bcftools concat -f ${OUT_DIR_RAW}/${trio}_file_list.txt -Oz -o ${OUT_DIR_RAW}/trio_${trio}.vcf.gz

  # --------------------------------------------
  # Step 3A: Prepare VCF for simulation
  # --------------------------------------------

  # Extract only EXOME variants
  bcftools view -i 'SNPSOURCE="EXOME"' -Oz \
    -o "${OUT_DIR_RAW}/trio_${trio}_exome.vcf.gz" "${OUT_DIR_RAW}/trio_${trio}.vcf.gz"

  # Keep only biallelic variants
  bcftools view -m2 -M2 -Oz \
    -o "${OUT_DIR_RAW}/trio_${trio}_exome_biallelic.vcf.gz" "${OUT_DIR_RAW}/trio_${trio}_exome.vcf.gz"

  # Apply region-based filtering
  bcftools view -T ^"$BED_REGIONS" -Oz \
    -o "${OUT_DIR_REGIONS}/filter_111_centro_seg_dup_trio_${trio}.vcf.gz" \
    "${OUT_DIR_RAW}/trio_${trio}_exome_biallelic.vcf.gz"

  # --------------------------------------------
  # Step 3B: Generate TSV file with SV coordinates for the proband
  # --------------------------------------------

  # Step 1: extract SVs from full trio VCF
  bcftools view -i 'VT="SV"' -Oz \
    -o "${OUT_DIR_RAW}/SV_trio_${trio}.vcf.gz" "${OUT_DIR_RAW}/trio_${trio}.vcf.gz"

  # Step 2: keep SVs where proband is not 0|0
  bcftools view -i 'GT!="0|0"' -s $proband -Oz \
    -o "${OUT_DIR_RAW}/${proband}_SV_filtered.vcf.gz" "${OUT_DIR_RAW}/SV_trio_${trio}.vcf.gz"

  # Step 3: generate raw TSV with SVLEN and END
  bcftools query -f '%CHROM\t%POS\t%INFO/SVLEN\n' \
    "${OUT_DIR_RAW}/${proband}_SV_filtered.vcf.gz" > \
    "${SV_COORD_DIR}/${trio}_SV_raw.tsv"

  # Step 4: extract final CHROM, START, END (computed from SVLEN)
  awk 'BEGIN{OFS="\t"}
    $3 != "." && $3 ~ /^-?[0-9]+$/ {
    svlen = ($3 < 0) ? -$3 : $3;
    print $1, $2, $2 + svlen;
    }' "${SV_COORD_DIR}/${trio}_SV_raw.tsv" > \
    "${SV_COORD_DIR}/${trio}_SV_filtered_not_INV.tsv"

  # Step 5: filter simulation-ready VCF by removing these SV regions
  bcftools view -T ^"${SV_COORD_DIR}/${trio}_SV_filtered_not_INV.tsv" -Oz \
    -o "${OUT_DIR_SV}/filter_111_centro_seg_dup_trio_${trio}_SV.vcf.gz" \
    "${OUT_DIR_REGIONS}/filter_111_centro_seg_dup_trio_${trio}.vcf.gz"

  echo "Done with $trio"

done

