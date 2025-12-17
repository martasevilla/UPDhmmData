#!/usr/bin/env Rscript

################################################################################
# Script: calculateEventsSfari.R
# Title: UPD event detection with UPDhmm on SFARI WGS trios
# Author: Marta Sevilla Porras
#
# Description:
# This script applies UPDhmm to trio-based WGS VCF files.
# The workflow includes:
#   - Sample validation and reordering using vcfCheck()
#   - VAF-based genotype correction using AD fields
#   - Quality filters (DP, GQ)
#   - Optional removal of centromeric, segmental duplication (SD)
#   - UPD event detection using calculateEvents()
#
# Notes:
#   Filtering of centro/SD regions can alternatively be performed
#   during VCF preprocessing using bcftools
#   (see 01b-createDataGenome.sh).
#
# Input arguments:
#   1) Chromosome
#   2) Input VCF directory
#   3) Output directory
#   4) Proband ID
#   5) Father ID
#   6) Mother ID
#   7) (Optional) TSV file with regions to exclude (SD / centromeres)
#
################################################################################

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(UPDhmm)
  library(future.apply)
  library(GenomicRanges)
  library(dplyr)
})

# -----------------------------
# Read arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    "Usage: calculateEventsSfari.R <chr> <input_vcf_dir> <output_dir> ",
    "<proband_id> <father_id> <mother_id> [exclude_regions.tsv]"
  )
}

chromosome  <- args[1]
input_dir   <- args[2]
output_dir  <- args[3]
proband_id  <- args[4]
father_id   <- args[5]
mother_id   <- args[6]
regions_tsv <- if (length(args) >= 7) args[7] else NA

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Parallelization
# -----------------------------
plan(multisession, workers = 8)

# -----------------------------
# Load exclusion regions (optional)
# -----------------------------
exclude_gr <- NULL

if (!is.na(regions_tsv) && file.exists(regions_tsv)) {
  message("Loading exclusion regions from: ", regions_tsv)
  reg_df <- read.table(regions_tsv, header = FALSE)
  colnames(reg_df)[1:3] <- c("chr", "start", "end")
  
  exclude_gr <- GRanges(
    seqnames = reg_df$chr,
    ranges   = IRanges(start = reg_df$start, end = reg_df$end)
  )
}

# -----------------------------
# VAF thresholds
# -----------------------------
vaf_low  <- 0.15
vaf_high <- 0.85

# -----------------------------
# Input VCF files
# -----------------------------
vcf_files <- list.files(input_dir, pattern = "\\.vcf.gz$", full.names = TRUE)

if (length(vcf_files) == 0) {
  stop("No VCF files found in input directory.")
}

# -----------------------------
# Main processing function
# -----------------------------
process_vcf <- function(vcf_file) {
  
  message("Processing: ", basename(vcf_file))
  
  vcf <- readVcf(vcf_file, genome = "hg38")
  
  # Validate and reorder trio using vcfCheck()
  vcf <- vcfCheck(
    vcf,
    proband = proband_id,
    father  = father_id,
    mother  = mother_id,
    run_qc  = FALSE
  )
  
  # Extract AD and compute VAF
  ad <- geno(vcf)$AD
  ns <- ncol(ad)
  nv <- nrow(ad)
  
  vaf <- matrix(NA, nrow = nv, ncol = ns)
  for (i in seq_len(nv)) {
    for (j in seq_len(ns)) {
      a <- ad[i, j][[1]]
      if (length(a) == 2 && sum(a) > 0) {
        vaf[i, j] <- a[2] / sum(a)
      }
    }
  }
  
  # Genotype correction
  gt <- geno(vcf)$GT
  gt[vaf < vaf_low]  <- "0/0"
  gt[vaf > vaf_high] <- "1/1"
  geno(vcf)$GT <- gt
  
  # Quality filters
  keep_dp <- apply(geno(vcf)$DP, 1, function(x) all(!is.na(x) & x >= 15))
  keep_gq <- apply(geno(vcf)$GQ, 1, function(x) all(!is.na(x) & x >= 80))
  keep_gt <- apply(geno(vcf)$GT, 1, function(x) any(x != "0/0"))
  
  vcf <- vcf[keep_dp & keep_gq & keep_gt]
  
  # Optional exclusion of SV / SD / centromeric regions
  if (!is.null(exclude_gr)) {
    vr <- rowRanges(vcf)
    hits <- findOverlaps(vr, exclude_gr, type = "within")
    keep <- setdiff(seq_along(vr), queryHits(hits))
    vcf <- vcf[keep]
  }
  
  # UPDhmm event detection
  events <- tryCatch(
    calculateEvents(vcf),
    error = function(e) {
      message("UPDhmm failed for ", basename(vcf_file), ": ", e$message)
      return(data.frame())
    }
  )
  
  return(events)
}

# -----------------------------
# Run analysis
# -----------------------------
results <- future_lapply(vcf_files, process_vcf)

out_file <- file.path(
  output_dir,
  paste0("UPDhmm_events_chr", chromosome, ".csv")
)

write.csv(bind_rows(results), out_file, row.names = FALSE)

message("UPDhmm analysis completed for chromosome ", chromosome)
