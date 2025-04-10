library(data.table)

# File paths and RSID column names
files <- list(
  HF      = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno1_EUR.tsv.gz",      col = "rsID"),
  HFrEF   = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno3_EUR.tsv.gz",      col = "rsID"),
  HFpEF   = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno4_EUR.tsv.gz",      col = "rsID"),
  DCM     = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/Jurgens_DCM_GWAS_META.tsv.gz",        col = "rsID"),
  HCM     = list(path = "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/GCST90435254.tsv",                    col = "rs_id"),
  LVEF    = list(path = "/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvef_filtered.tsv.gz",  col = "SNP"),
  LVESVi  = list(path = "/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvesvi_filtered.tsv.gz", col = "SNP"),
  LVEDVi  = list(path = "/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvedvi_filtered.tsv.gz", col = "SNP"),
  LVMi    = list(path = "/medpop/afib/skhurshid/lvm_gwas/gwas/v19_seg/v19_seg_lvmi_adjusted_bolt.imputed.filtered.tsv", col = "SNP")
)

# Initialize the intersection set
rsids_in_all <- NULL

# Iterate and intersect RSIDs
for (i in seq_along(files)) {
  file_info <- files[[i]]
  message("Processing: ", file_info$path)
  dt <- fread(file_info$path, select = file_info$col, showProgress = FALSE)
  rsids <- dt[[1]]
  
  if (is.null(rsids_in_all)) {
    rsids_in_all <- rsids
  } else {
    rsids_in_all <- intersect(rsids_in_all, rsids)
  }
}

# Save result
fwrite(data.table(rsid = rsids_in_all), "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/rsids_in_all_files.csv")
