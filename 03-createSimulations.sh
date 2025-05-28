#!/bin/bash

################################################################################
# Script: 03-createSimulations.sh
# Title: Simulation of UPD Events for Benchmarking (Exome or Genome)
# Author: Marta Sevilla Porras
# Date: 26/05/2025
#
# Description:
# Simulates exactly 100 UPD events in total across real trio VCFs.
# Each event applies a UPD type (iso/het × maternal/paternal) to a region of
# the proband’s genome. All 20 combinations of 4 UPD types × 5 region sizes
# are represented, each repeated 5 times across the available trios.
#
# Inputs:
# - config/trios.tsv                 → Trio list: TrioID, Proband, Father, Mother
# - regions/filtered_final_<size>.txt → Region lists per size
# - data/01_centro_seg_dup_SV/      → Filtered VCFs per trio (SVs)
# - scripts/processRegion.R         → R script that modifies the proband genotype
#
# Outputs:
# - Simulated VCFs in:
#   results_simulation/updhmm/<size>/
#   results_simulation/updio/<size>/
#   results_simulation/updbatch/<size>/
#
# Usage:
#   sbatch --array=1-100 03-createSimulations.sh
################################################################################

#SBATCH --array=1-100

# --------------------------------------------
# Configuration
# --------------------------------------------


TRIO_FILE="config/trios.tsv"
SCRIPT_R="scripts/processRegion.R"
REGION_SIZES=("1mb" "2mb" "5mb" "10mb" "20mb")
UPD_TYPES=("het_fat" "het_mat" "iso_fat" "iso_mat")

NUM_SIMULATIONS=100
NUM_COMBINATIONS=20  # 4 UPD types × 5 region sizes

module load R/4.1.2-foss-2021b
module load BCFtools/1.16-GCC-11.2.0

# --------------------------------------------
# Load trios
# --------------------------------------------
readarray -t trios < <(tail -n +1 "$TRIO_FILE")
NUM_TRIOS=${#trios[@]}

if (( NUM_TRIOS == 0 )); then
  echo "No trios found in $TRIO_FILE"
  exit 1
fi

# --------------------------------------------
# Determine this task's UPD type, region size, and trio
# --------------------------------------------
task_id=$SLURM_ARRAY_TASK_ID
combo_index=$(( (task_id - 1) % NUM_COMBINATIONS ))    # 0 to 19
repeat_index=$(( (task_id - 1) / NUM_COMBINATIONS ))   # 0 to 4

upd_type="${UPD_TYPES[$(( combo_index / 5 ))]}"
region_size="${REGION_SIZES[$(( combo_index % 5 ))]}"
trio_index=$(( task_id % NUM_TRIOS ))
[[ $trio_index -eq 0 ]] && trio_index=$NUM_TRIOS
trio_index=$(( trio_index - 1 ))

# Trio info
trio_line="${trios[$trio_index]}"
IFS=$'\t' read -r TRIO_ID PROBAND FATHER MOTHER <<< "$trio_line"
samples=("$PROBAND" "$FATHER" "$MOTHER")

# Region file and region line
REGION_FILE="regions/filtered_final_${region_size}.txt"
region=$(sed -n "${task_id}p" "$REGION_FILE")
IFS=$'\t' read -r chr pos1 pos2 _ <<< "$region"
region_tag="${chr}_${pos1}"


vcf_input="data/01_centro_seg_dup_SV/filter_111_centro_seg_dup_trio_${TRIO_ID}_SV.vcf.gz"

tabix "$vcf_input"

# --------------------------------------------
# Extract and simulate region
# --------------------------------------------
bcftools view -r "$chr:$pos1-$pos2" -Oz -o "${region_tag}_selected_region.vcf.gz" "$vcf_input"
bcftools view -t ^"$chr:$pos1-$pos2" -Oz -o "${region_tag}_excluded_region.vcf.gz" "$vcf_input"
tabix "${region_tag}_selected_region.vcf.gz"
tabix "${region_tag}_excluded_region.vcf.gz"

# Run simulation
Rscript "$SCRIPT_R" "$upd_type" "$chr" "$pos1"
bgzip "${region_tag}_edited.vcf"
tabix "${region_tag}_edited.vcf.gz"

# Merge VCF
merged_unsorted="tmp_${region_tag}_merged.vcf.gz"
merged_sorted="${TRIO_ID}_${upd_type}_${region_tag}_merge_sort.vcf.gz"

bcftools concat -a "${region_tag}_edited.vcf.gz" "${region_tag}_excluded_region.vcf.gz" -Oz -o "$merged_unsorted"
bcftools sort "$merged_unsorted" -Oz -o "$merged_sorted"
tabix "$merged_sorted"

# --------------------------------------------
# Output directories
# --------------------------------------------
dir_updhmm="results_simulation/updhmm/${region_size}/"
dir_updio="results_simulation/updio/${region_size}/"
dir_altaf="results_simulation/updbatch/${region_size}/"

mkdir -p "$dir_updhmm" "$dir_updio" "$dir_altaf"
cp "$merged_sorted" "${dir_updhmm}/${merged_sorted}"
cp "$merged_sorted" "${dir_updio}/${merged_sorted}"

# --------------------------------------------
# Generate individual VCFs (AltAFplotter)
# --------------------------------------------
roles=("proband" "father" "mother")

for i in {0..2}; do
  sample="${samples[$i]}"
  role="${roles[$i]}"
  indiv_base="${TRIO_ID}_${upd_type}_${region_tag}_${role}"

  bcftools view -s "$sample" "$merged_sorted" -Oz -o "${dir_altaf}/${indiv_base}_prov.vcf.gz"
  tabix "${dir_altaf}/${indiv_base}_prov.vcf.gz"

  bcftools view -e 'FORMAT/GT=="0|0"' "${dir_altaf}/${indiv_base}_prov.vcf.gz" -Oz -o "${dir_altaf}/${indiv_base}_merge_sort.vcf.gz"
  tabix "${dir_altaf}/${indiv_base}_merge_sort.vcf.gz"
  rm "${dir_altaf}/${indiv_base}_prov.vcf.gz"*
done

echo "Simulation complete: $TRIO_ID | $upd_type | $region_size | $chr:$pos1"