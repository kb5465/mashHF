################################################################################
# Title: Subset GWAS Files to eQTL SNPs
# Description:
#   Filters multiple GWAS summary statistic files to retain only SNPs that are
#   significant eQTLs in GTEx v8 (as precompiled in `mashHF_eqtl_031825.csv`).
#   Outputs harmonized GWAS subsets for use in downstream colocalization or 
#   gene-trait association testing.
#
# Author: Kiran Biddinger
# Date: April 2025
################################################################################

library(data.table)

# ----------------------------
# Load eQTL SNP list
# ----------------------------
eqtls <- fread("/medpop/afib/kbiddinger/Projects/mashHF/data/mashHF_eqtl_031825.csv")
snp_list <- eqtls$snpId

# ----------------------------
# Output directory
# ----------------------------
output_dir <- "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/"

# ----------------------------
# Helper function
# ----------------------------
process_gwas <- function(file, snp_col, a1_col, a2_col, beta_col, se_col, outname) {
  gwas <- fread(file)
  gwas <- subset(gwas, get(snp_col) %in% snp_list)
  gwas <- gwas[, .SD, .SDcols = c(snp_col, a1_col, a2_col, beta_col, se_col)]
  setnames(gwas, c("rsid", "ea", "nea", "beta", "se"))
  fwrite(gwas, paste0(output_dir, outname, ".csv"))
}

# ----------------------------
# Run for each trait
# ----------------------------
process_gwas("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno1_EUR.tsv.gz",
             "rsID", "A1", "A2", "A1_beta", "se", "HF")

process_gwas("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno3_EUR.tsv.gz",
             "rsID", "A1", "A2", "A1_beta", "se", "HFrEF")

process_gwas("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno4_EUR.tsv.gz",
             "rsID", "A1", "A2", "A1_beta", "se", "HFpEF")

process_gwas("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/Jurgens_DCM_GWAS_META.tsv.gz",
             "rsID", "EA", "NEA", "BETA", "SE", "DCM")

process_gwas("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/GCST90435254.tsv",
             "rs_id", "effect_allele", "other_allele", "beta", "standard_error", "HCM")

process_gwas("/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvef_filtered.tsv.gz",
             "SNP", "ALLELE1", "ALLELE0", "BETA", "SE", "LVEF")

process_gwas("/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvesvi_filtered.tsv.gz",
             "SNP", "ALLELE1", "ALLELE0", "BETA", "SE", "LVESVi")

process_gwas("/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvedvi_filtered.tsv.gz",
             "SNP", "ALLELE1", "ALLELE0", "BETA", "SE", "LVEDVi")

process_gwas("/medpop/afib/skhurshid/lvm_gwas/gwas/v19_seg/v19_seg_lvmi_adjusted_bolt.imputed.filtered.tsv",
             "SNP", "ALLELE1", "ALLELE0", "BETA", "SE", "LVMi")