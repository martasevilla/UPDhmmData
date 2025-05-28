#!/bin/bash

################################################################################
# Script: 00-processSFARIrawVCFs.sh
# Title: Process raw WGS SFARI cohort files for UPDhmm event detection
# Author: Marta Sevilla Porras
# Date: 29/05/2025
#
# Description:
# This script illustrates the workflow for whole-genome
# sequencing data from the SFARI cohort. The process includes:
#  - Extraction of variants (SNVs and indels)
#  - Removal of unnecessary annotations and formats
#  - Extraction and combination of structural variant (SV) regions per trio
#  - Execution of an R script to calculate UPD events using UPDhmm with
#    VAF-based genotype corrections and quality filters
#
# Requirements:
#  - bcftools
#  - R packages: UPDhmm, VariantAnnotation, future.apply, GenomicRanges, dplyr
################################################################################

# -----------------------------
# Step 1: Extract SNVs and indels from the raw VCF file
# -----------------------------
bcftools view -v snps,indels \
  -Ob -o 01_filter/filter_remove_chr15_pre.bcf \
  00_raw/CCDG_9000JG_B01_GRM_WGS_2019-03-21_chr15.recalibrated_variants.masked.vcf.gz

# -----------------------------
# Step 2: Retain only relevant FORMAT fields; remove all INFO fields
# -----------------------------
bcftools annotate \
  -x INFO,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/AB,^FORMAT/FT,^FORMAT/GQ,^FORMAT/PL \
  -Ob -o 01_filter/filter_remove_chr15_batch6.bcf \
  01_filter/filter_remove_chr15_pre.bcf

# -----------------------------
# Step 3: Extract a specific trio (example: SSC03417 trio)
# Filters:
#   - Keep only positions where not all three samples are 0/0
#   - Keep only biallelic variants
#   - Remove missing genotypes
# -----------------------------
bcftools view \
  -s SSC03417,SSC03430,SSC03423 \
  01_filter/filter_remove_chr15_batch6.bcf \
  -i 'COUNT(FORMAT/GT="0/0") != 3' \
  -m2 -M2 -e 'GT="./."' \
  -Oz -o 02_vcf_files_trios/chr15/chr15_filter_111_SSC03417_remove_SV.vcf.gz

# -----------------------------
# Step 4: Generate SV region files per individual and combine into trio-level TSVs
# The SFARI project is divided into multiple projects and sequencing batches.
# For each individual, structural variant (SV) coordinates have been extracted separately.
# Then, for every trio, the SV files of the proband, father, and mother are merged into a single TSV file.
# These merged files are stored in:
#   SV/trios_SV/<proband_ID>_trio.tsv
# -----------------------------


# -----------------------------
# Step 5: Calculate UPD events using UPDhmm with VAF-based corrections
# This step will apply the `calculateEventsSfari.R` script to all VCFs
# in the specified chromosome directory (example: chr15).
# -----------------------------
Rscript calculateEventsSfari.R 15
