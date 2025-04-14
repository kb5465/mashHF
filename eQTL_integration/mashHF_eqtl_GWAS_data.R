library(data.table)

# Read eQTL data
eqtls <- fread("/medpop/afib/kbiddinger/Projects/mashHF/data/mashHF_eqtl_031825.csv")

# HF
gwas <- fread("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno1_EUR.tsv.gz")
gwas <- subset(gwas, rsID %in% eqtls$snpId)
gwas <- gwas[,c("rsID", "A1", "A2", "A1_beta", "se")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/HF.csv")

# HFrEF
gwas <- fread("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno3_EUR.tsv.gz")
gwas <- subset(gwas, rsID %in% eqtls$snpId)
gwas <- gwas[,c("rsID", "A1", "A2", "A1_beta", "se")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/HFrEF.csv")

# HFpEF
gwas <- fread("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/FORMAT-METAL_Pheno4_EUR.tsv.gz")
gwas <- subset(gwas, rsID %in% eqtls$snpId)
gwas <- gwas[,c("rsID", "A1", "A2", "A1_beta", "se")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/HFpEF.csv")

# DCM GWAS
gwas <- fread("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/Jurgens_DCM_GWAS_META.tsv.gz")
gwas <- subset(gwas, rsID %in% eqtls$snpId)
gwas <- gwas[,c("rsID", "EA", "NEA", "BETA", "SE")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/DCM.csv")

# HCM GWAS
gwas <- fread("/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS/GCST90435254.tsv")
gwas <- subset(gwas, rs_id %in% eqtls$snpId)
gwas <- gwas[,c("rs_id", "effect_allele", "other_allele", "beta", "standard_error")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/HCM.csv")

# LVEF GWAS
gwas <- fread("/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvef_filtered.tsv.gz")
gwas <- subset(gwas, SNP %in% eqtls$snpId)
gwas <- gwas[,c("SNP", "ALLELE1", "ALLELE0", "BETA", "SE")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/LVEF.csv")

# LVESVi GWAS
gwas <- fread("/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvesvi_filtered.tsv.gz")
gwas <- subset(gwas, SNP %in% eqtls$snpId)
gwas <- gwas[,c("SNP", "ALLELE1", "ALLELE0", "BETA", "SE")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/LVESVi.csv")

# LVEDVi GWAS
gwas <- fread("/medpop/esp2/karagam/UK_Biobank/phenos/HF_500K/Cardiac_MRI/cMRI_36K_v68/lvedvi_filtered.tsv.gz")
gwas <- subset(gwas, SNP %in% eqtls$snpId)
gwas <- gwas[,c("SNP", "ALLELE1", "ALLELE0", "BETA", "SE")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/LVEDVi.csv")

# LVMi GWAS
gwas <- fread("/medpop/afib/skhurshid/lvm_gwas/gwas/v19_seg/v19_seg_lvmi_adjusted_bolt.imputed.filtered.tsv")
gwas <- subset(gwas, SNP %in% eqtls$snpId)
gwas <- gwas[,c("SNP", "ALLELE1", "ALLELE0", "BETA", "SE")]
colnames(gwas) <- c("rsid", "ea", "nea", "beta", "se")
write.csv(gwas, "/medpop/afib/kbiddinger/Projects/mashHF/data/GWAS_subset/LVMi.csv")
