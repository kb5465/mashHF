#!/usr/bin/env Rscript

#===============================================================================
# Script: intersect_rsids_across_gwas.R
# Author: Kiran Biddinger
# Date: 2025-04-14
#
# Description:
#   This script identifies the set of SNP rsIDs that are shared across all GWAS 
#   summary statistics used in the mashHF project. It performs an intersection 
#   over multiple GWAS datasets and saves the common set of SNPs to a file.
#
# Usage:
#   Run this script with Rscript:
#     Rscript intersect_rsids_across_gwas.R
#
# Output:
#   A file listing SNP rsIDs present in all input GWAS datasets:
#     /medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/rsids_in_all_files.csv
#
# Dependencies:
#   - data.table
#
# Notes:
#   - GWAS files must contain rsID columns as specified below.
#   - Adjust file paths or column names as needed if sources change.
#===============================================================================

library(data.table)

# Define GWAS file paths and RSID column names
gwas_files <- list(
  HF      = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno1_EUR.tsv.gz",       col = "rsID"),
  HFrEF   = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno3_EUR.tsv.gz",       col = "rsID"),
  HFpEF   = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno4_EUR.tsv.gz",       col = "rsID"),
  DCM     = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/Jurgens_DCM_GWAS_META.tsv.gz",         col = "rsID"),
  HCM     = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/GCST90435254.tsv",                     col = "rs_id"),
  LVEF    = list(path = "/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvef_filtered.tsv.gz",   col = "SNP"),
  LVESVi  = list(path = "/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvesvi_filtered.tsv.gz", col = "SNP"),
  LVEDVi  = list(path = "/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvedvi_filtered.tsv.gz", col = "SNP"),
  LVMi    = list(path = "/medpop/afib/skhurshid/lvm_gwas/gwas/v19_seg/v19_seg_lvmi_adjusted_bolt.imputed.filtered.tsv",  col = "SNP")
)

# Initialize set of common RSIDs
rsids_in_all <- NULL

# Iterate through files to get intersection
for (trait in names(gwas_files)) {
  file_info <- gwas_files[[trait]]
  message("Processing ", trait, " from: ", file_info$path)
  
  dt <- fread(file_info$path, select = file_info$col, showProgress = FALSE)
  rsids <- dt[[1]]
  
  # Update intersection
  rsids_in_all <- if (is.null(rsids_in_all)) rsids else intersect(rsids_in_all, rsids)
}

# Output path
output_path <- "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/rsids_in_all_files.csv"

# Write intersected RSIDs to CSV
fwrite(data.table(rsid = rsids_in_all), output_path)

message("Saved intersected RSIDs to: ", output_path)
