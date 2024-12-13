#!/bin/bash

# Define phenotype list
PHENO_LIST=(
    HF
    DCM
    HCM
    NTproBNP
    Troponin
    LVMi
    LVESVi
    LVEDVi
    LVEF
    LVConc
    meanWT
    strain_rad
    strain_long
    strain_circ
)

# Iterate over PHENO_LIST
for PHENO in "${PHENO_LIST[@]}"; do
    # Directory containing the files
    PHENO_DIR="/mnt/project/kbiddinger/projects/mashHF/results/step2/rawassociation/${PHENO}"

    # Output directory and file
    OUTPUT_DIR="exome-seq:/kbiddinger/projects/mashHF/results/prior/"
    OUTPUT_FILE="${PHENO}_canonical_LOF_0.01_merged.tsv"

    # Create the output directory if it doesn't exist
    dx mkdir -p "${OUTPUT_DIR}"

    # Write file header
    echo -e "ID\tBETA_${PHENO}\tSE_${PHENO}" > "${OUTPUT_FILE}"

    # Loop through all files in the PHENO directory
    for file in "${PHENO_DIR}"/*/*; do
        # Process the file
        awk '
        BEGIN { OFS="\t" }
        # Match rows where ID contains "canonical.LOF.0.01" and TEST is exactly "ADD"
        $3 ~ /canonical\.LOF\.0\.01/ && $8 == "ADD" {
            # Extract ENSG ID by splitting the $3 field
            split($3, parts, "_")
            ensID = parts[1]  # Assume the ENSG ID is the first part
            if (ensID ~ /^ENSG[0-9]+$/) {
                print ensID, $9, $10
            }
        }' "$file" >> "${OUTPUT_FILE}"
    done

    # Upload the output file to the specified directory
    dx upload "${OUTPUT_FILE}" --destination "${OUTPUT_DIR}"
done