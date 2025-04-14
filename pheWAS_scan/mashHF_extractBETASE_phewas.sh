#!/bin/bash

# NOTE: Last time outputted to exome-seq:/kbiddinger/projects/mashHF/results/, even though
# `prior` directory was made

# Define phenotype list
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

# Iterate over PHENO_LIST
for PHENO in "${PHENO_LIST[@]}"; do
    # Directory containing the files
    PHENO_DIR="/mnt/project/kbiddinger/projects/mashHF/results/step2phewas/rawassociation/${PHENO}"

    # Output directory and file
    OUTPUT_DIR="exome-seq:/kbiddinger/projects/mashHF/results/st_phewas/"
    OUTPUT_FILE="${PHENO}_canonical_LOF_0.01_merged.tsv"

    # Create the output directory if it doesn't exist
    dx mkdir -p "${OUTPUT_DIR}"

    # Write file header
    echo -e "ID\tBETA_${PHENO}\tSE_${PHENO}\tN_CARRIERS\tN_TOTAL" > "${OUTPUT_FILE}"

    # Loop through all files in the PHENO directory
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
            }' "$file" >> "$OUTPUT_FILE"
        fi
    done

    # Upload the output file to the specified directory
    dx rm "${OUTPUT_DIR}"/"${OUTPUT_FILE}"
    dx upload "${OUTPUT_FILE}" --destination "${OUTPUT_DIR}"
done