#!/bin/bash

################################################################################
# Script: 01-createData_genome.sh
# Title: Preprocessing Genome VCFs for Trio-Based UPD Simulations
# Author: Marta Sevilla Porras
# Date: 27/05/2025
#
# Description:
# This script preprocesses high-coverage genome VCFs from the 1000 Genomes Project.
# It:
#   - Extracts 10 trios from per-chromosome VCFs
#   - Concatenates them into full-genome VCFs
#   - Applies region-based filtering (centromeres, segdups, HLA/KIR)
#   - Removes fully homozygous ref positions
#   - Optionally filters SV regions (non-inversions)
#   - Produces final simulation-ready VCFs and SV exclusion BEDs
#
# Output VCFs:
#   - data/01_centro_seg_dup/         (region-only filtered)
#   - data/01_centro_seg_dup_SV/      (region + SV filtered)
#   - SVs_coordinates/<trio>_SV_filtered_not_INV.tsv
#
# Requirements:
#   - bcftools ≥ 1.16
# 
################################################################################

# -----------------------------
# Configuration
# -----------------------------
TRIOS=("HG00405  HG00429  HG00552  HG01674  HG01881  HG01981  HG03008  HG03296  HG03605  NA12386  NA18503")
FATHER="father" #should be filled with correspondence father ID
MOTHER="mother" #should be filled with correspondence mother ID
CHROMOSOMES=$(seq 1 22)
RAW_BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
BED_REGIONS="centro_seg_dup_coordinates/centro_seg_dup.tsv"
SV_BASE="SVs_coordinates"
RAW_DIR="data/00_raw"
REGION_OUT="data/01_centro_seg_dup"
SV_OUT="data/01_centro_seg_dup_SV"

mkdir -p "$RAW_DIR" "$REGION_OUT" "$SV_OUT" "$SV_BASE"

# -----------------------------
# Process each trio
# -----------------------------
for TRIO_ID in "${TRIOS[@]}"; do
  echo "Processing trio: $TRIO_ID"

  # -----------------------------
  # Step 1: Download and extract trio VCFs per chromosome
  # -----------------------------
  for chr in $CHROMOSOMES; do
    bcftools view \
      -s ${TRIO_ID},${FATHER},${MOTHER} \
      --force-samples \
      -Oz -o "${RAW_DIR}/${TRIO_ID}_vcf_chr${chr}.vcf.gz" \
      "${RAW_BASE_URL}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

    bcftools annotate --remove INFO,ID,^FORMAT/GT \
      -Oz -o "${RAW_DIR}/${TRIO_ID}_edited_chr${chr}.vcf.gz" \
      "${RAW_DIR}/${TRIO_ID}_vcf_chr${chr}.vcf.gz"
  done

  # -----------------------------
  # Step 2: Concatenate chromosomes
  # -----------------------------
  ls "${RAW_DIR}/${TRIO_ID}_edited_chr"*.vcf.gz > "${RAW_DIR}/${TRIO_ID}_file_list.txt"
  bcftools concat -f "${RAW_DIR}/${TRIO_ID}_file_list.txt" -Oz -o "${RAW_DIR}/concat_${TRIO_ID}.vcf.gz"

  
  
  # --------------------------------------------
  # Step 3A: Prepare VCF for simulation
  # --------------------------------------------

# Apply region-based filtering
  bcftools view -v snps,indels -T ^"$BED_REGIONS" -Oz \
    -o "${OUT_DIR_REGIONS}/filter_centro_seg_dup_trio_concat_${TRIO_ID}.vcf.gz" \
    "${RAW_DIR}/concat_${TRIO_ID}.vcf.gz"

  # Remove positions where all samples are homozygous reference
  bcftools view -i 'COUNT(GT="0|0")!=3' -Oz \
    -o "${OUT_DIR_REGIONS}/filter_111_centro_seg_dup_trio_${TRIO_ID}.vcf.gz" \
    "${OUT_DIR_REGIONS}/filter_centro_seg_dup_trio_concat_${TRIO_ID}.vcf.gz"



  # --------------------------------------------
  # Step 3B: Generate TSV file with SV coordinates for the proband
  # --------------------------------------------

# Step 1: extract SVs where proband is not 0|0
  bcftools view -v other -i 'GT!="0|0"'  -Oz \
    -o "${RAW_DIR}/${proband}_SV_filtered.vcf.gz" "${RAW_DIR}/concat_${TRIO_ID}.vcf.gz"

  # Step 2: remove inversions (SVTYPE="INV")
  bcftools view -e 'SVTYPE="INV"' -Oz \
    -o "${RAW_DIR}/${proband}_SV_filtered_not_INV.vcf.gz" \
    "${RAW_DIR}/${proband}_SV_filtered.vcf.gz"

  # Step 3: generate raw TSV with SVLEN and END
  bcftools query -f '%CHROM\t%POS\t%INFO/SVLEN\n' \
    "${RAW_DIR}/${proband}_SV_filtered_not_INV.vcf.gz" > \
    "${SV_COORD_DIR}/${TRIO_ID}_SV_raw.tsv"

  # Step 4: extract final CHROM, START, END (computed from SVLEN)
  awk 'BEGIN{OFS="\t"} \
    $3 != "." && $3 ~ /^-?[0-9]+$/ {
      svlen = ($3 < 0) ? -$3 : $3;
      print $1, $2, $2 + svlen;
    }' "${SV_COORD_DIR}/${TRIO_ID}_SV_raw.tsv" > \
    "${SV_COORD_DIR}/${TRIO_ID}_SV_filtered_not_INV.tsv"

  # Step 5: filter simulation-ready VCF by removing these SV regions
  bcftools view -T ^"${SV_COORD_DIR}/${TRIO_ID}_SV_filtered_not_INV.tsv" -Oz \
    -o "${OUT_DIR_SV}/filter_111_centro_seg_dup_trio_${TRIO_ID}_SV.vcf.gz" \
    "${OUT_DIR_REGIONS}/filter_111_centro_seg_dup_trio_${TRIO_ID}.vcf.gz"

  echo "Done with $TRIO_ID"

done
