#!/bin/bash

###############################################
# mashHF: Extract and Merge Canonical LOF 0.01 Results
# Author: Kiran Biddinger
# Description:
#   - For each phenotype, merges all gene-level association
#     results matching `canonical.LOF.0.01` from REGENIE step 2
#   - Calculates carrier count = allele frequency × N
#   - Writes standardized output file and uploads to DNAnexus
###############################################

# ----------------------------
# 1. Phenotype List
# ----------------------------
PHENO_LIST=(
  HF DCM HCM
  NTproBNP Troponin
  LVMi LVESVi LVEDVi LVEF LVConc
  meanWT strain_rad strain_long strain_circ
)

# ----------------------------
# 2. Loop Through Each Phenotype
# ----------------------------
for PHENO in "${PHENO_LIST[@]}"; do

  # Define input and output locations
  PHENO_DIR="/mnt/project/kbiddinger/projects/mashHF/results/step2/rawassociation/${PHENO}"
  OUTPUT_DIR="exome-seq:/kbiddinger/projects/mashHF/results/st/"
  OUTPUT_FILE="${PHENO}_canonical_LOF_0.01_merged.tsv"

  # Create destination directory (if not already present)
  dx mkdir -p "${OUTPUT_DIR}"

  # Remove previous version of file, if it exists
  dx rm "${OUTPUT_DIR}/${OUTPUT_FILE}"

  # Initialize output file with header
  echo -e "ID\tBETA_${PHENO}\tSE_${PHENO}\tN_CARRIERS\tN_TOTAL" > "${OUTPUT_FILE}"

  # ----------------------------
  # 3. Process Files for This Phenotype
  # ----------------------------
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

  # ----------------------------
  # 4. Upload Merged Result File
  # ----------------------------
  dx upload "${OUTPUT_FILE}" --destination "${OUTPUT_DIR}"

done

# ----------------------------
# 5. Clean Workspace
# ----------------------------
rm -rf /opt/notebooks/*