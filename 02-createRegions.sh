#!/bin/bash

################################################################################
# Script: 02-createRegions.sh
# Title: Simulation Region Selector for UPD Simulations (Exome or Genome)
# Author: Marta Sevilla Porras
# Date: 26/05/2025
#
# Description:
# This script identifies suitable genomic intervals for simulation based on:
#   - A fixed interval size (e.g. 1Mb, 2Mb, 5Mb, 10Mb, 20Mb)
#   - Real exome data with sufficient variant density
#
# Steps:
#   1. Generate random genomic intervals
#   2. Intersect each VCF with those intervals
#   3. Filter intervals with at least 10 variants per sample
#   4. Select 100 valid intervals randomly
#   5. Format results into BED-like files for downstream use
#
# Requirements:
#   - BEDTools (≥ 2.29)
#   - bcftools (≥ 1.19)
#
# Input:
#   - genome.txt: Chromosome sizes (BEDTools format)
#   - vcf_files.txt: List of pre-filtered .vcf.gz files (indexed with tabix)
#
# Output:
#   - BED file with 100 selected regions per size (e.g., regions_1000000_exome.bed)
################################################################################

# -----------------------------
# Parameters (Adjustable)
# -----------------------------

INTERVAL_SIZE=1000000                # 1Mb. Change as needed for 2Mb, 5Mb, etc.
NUM_INTERVALS=10000                  # Initial number of intervals to test
GENOME_FILE="genome.txt"            # Genome chromosome sizes
VCF_LIST="vcf_files.txt"            # List of input VCFs (tabix-indexed)

# Derived filenames
PREFIX="regions_${INTERVAL_SIZE}"
INTERVALS_FILE="${PREFIX}.bed"
OUTPUT_DIR="intersect_${INTERVAL_SIZE}"
MERGED_COUNTS="${PREFIX}_merged.txt"
FILTERED="${PREFIX}_filtered.txt"
FORMATTED_BED="${PREFIX}_exome.bed"

# -----------------------------
# Step 1: Generate Random Intervals
# -----------------------------
echo "Generating ${INTERVAL_SIZE}-bp intervals..."
bedtools random -l $INTERVAL_SIZE -n $NUM_INTERVALS -g $GENOME_FILE > $INTERVALS_FILE

# -----------------------------
# Step 2: Intersect Each VCF
# -----------------------------
echo "Intersecting VCFs with intervals..."
mkdir -p $OUTPUT_DIR

while IFS= read -r vcf; do
    output_file="$OUTPUT_DIR/$(basename "$vcf")_intersect.bed"
    bedtools intersect -a "$INTERVALS_FILE" -b "$vcf" -c > "$output_file"
done < "$VCF_LIST"

# -----------------------------
# Step 3: Merge Intersections
# -----------------------------
echo "Merging counts..."
TMP="$OUTPUT_DIR/merged_tmp.txt"
first_file=$(head -n 1 "$VCF_LIST")
cut -f 1-3 "$OUTPUT_DIR/$(basename "$first_file")_intersect.bed" > "$TMP"

while IFS= read -r vcf; do
    col=$(cut -f 4 "$OUTPUT_DIR/$(basename "$vcf")_intersect.bed")
    paste "$TMP" <(echo "$col") > "${TMP}_new"
    mv "${TMP}_new" "$TMP"
done < "$VCF_LIST"

mv "$TMP" "$MERGED_COUNTS"

# -----------------------------
# Step 4: Filter by Variant Count
# -----------------------------
echo "Filtering intervals with at least 10 variants per sample..."
awk '
BEGIN { OFS="\t" }
{
  keep = 1;
  for (i = 4; i <= NF; i++) {
    if ($i < 10) {
      keep = 0;
      break;
    }
  }
  if (keep) print;
}
' "$MERGED_COUNTS" > "$FILTERED"

# -----------------------------
# Step 5: Select 100 Random Regions
# -----------------------------
echo "Selecting 100 random valid intervals..."
shuf -n 100 "$FILTERED" > "${FILTERED}.tmp" && mv "${FILTERED}.tmp" "$FILTERED"

# -----------------------------
# Step 6: Format as BED-like Table
# -----------------------------
echo "Formatting output..."
awk -v SIZE=$INTERVAL_SIZE '
BEGIN { OFS="\t"; print "seqnames", "start", "end", "width", "strand", "pos", "size" }
{
  seqnames = $1;
  start = $2;
  end = $3;
  width = end - start + 1;
  strand = (NR % 2 == 0) ? "+" : "-";
  pos = NR;
  print seqnames, start, end, width, strand, pos, SIZE;
}
' "$FILTERED" > "$FORMATTED_BED"

echo "Output written to: $FORMATTED_BED"

################################################################################
# To repeat for 2Mb, 5Mb, 10Mb, 20Mb:
#   ➤ Change INTERVAL_SIZE and re-run the script
#   ➤ Output will be named accordingly
################################################################################