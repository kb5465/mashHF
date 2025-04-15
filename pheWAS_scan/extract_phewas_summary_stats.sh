#!/bin/bash

################################################################################
# Script: extract_phewas_summary_stats.sh
# Author: Kiran Biddinger
# Project: mashHF
#
# Description:
#   Extracts gene-level BETA, SE, and carrier counts from REGENIE Step 2
#   results (rawassociation) for canonical LOF (0.01 max AAF) masks
#   across binary and quantitative PHEWAS traits.
#
# Usage:
#   bash extract_phewas_summary_stats.sh
#
# Inputs:
#   - REGENIE Step 2 output files (all chromosomes, chunks) in:
#       /mnt/project/kbiddinger/projects/mashHF/results/step2phewas/rawassociation/${PHENO}
#
# Outputs:
#   - Tab-delimited summary file for each phenotype:
#       ${PHENO}_canonical_LOF_0.01_merged.tsv
#     uploaded to:
#       exome-seq:/kbiddinger/projects/mashHF/results/st_phewas/
################################################################################

#------------------------------------------------------------------------------
# Define PHEWAS traits
#------------------------------------------------------------------------------

PHENO_LIST=(
    CAD
    T2D
    AF
    CKD
    IS
    LDL_C
    BMI
    SBP
    DBP
)

#------------------------------------------------------------------------------
# Process each phenotype
#------------------------------------------------------------------------------

for PHENO in "${PHENO_LIST[@]}"; do
    PHENO_DIR="/mnt/project/kbiddinger/projects/mashHF/results/step2phewas/rawassociation/${PHENO}"
    OUTPUT_DIR="exome-seq:/kbiddinger/projects/mashHF/results/st_phewas/"
    OUTPUT_FILE="${PHENO}_canonical_LOF_0.01_merged.tsv"

    # Create remote output directory (if needed)
    dx mkdir -p "${OUTPUT_DIR}"

    # Write header
    echo -e "ID\tBETA_${PHENO}\tSE_${PHENO}\tN_CARRIERS\tN_TOTAL" > "${OUTPUT_FILE}"

    # Loop through all .regenie files in PHENO_DIR
    for file in "${PHENO_DIR}"/*/*; do
        if [[ -f "$file" ]]; then
            awk '
            BEGIN { OFS="\t" }
            $3 ~ /canonical\.LOF\.0\.01/ && $8 == "ADD" {
                split($3, parts, "_")
                ensID = parts[1]
                a1freq = $6
                n = $7
                if (ensID ~ /^ENSG[0-9]+$/ && a1freq != "NA" && n != "NA") {
                    carriers = a1freq * n
                    print ensID, $9, $10, carriers, n
                }
            }' "$file" >> "${OUTPUT_FILE}"
        fi
    done

    # Remove previous version and upload new summary file
    dx rm "${OUTPUT_DIR}/${OUTPUT_FILE}"
    dx upload "${OUTPUT_FILE}" --destination "${OUTPUT_DIR}"
done