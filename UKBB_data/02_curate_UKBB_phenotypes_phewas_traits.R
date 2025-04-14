###############################################
# mashHF: PheWAS-Style UKBB Phenotype Curation
# Author: Kiran Biddinger
# Description: Prepares a clean dataset with 
# cardiometabolic biomarkers and binary traits 
# for downstream PheWAS-style analyses.
###############################################

# =======================
# 1. Load Required Libraries
# =======================
library(data.table)
library(dplyr)

# =======================
# 2. Load Bulk UKBB Data and Select Variables
# =======================
ukbb <- fread("/mnt/project/sjj/short_projects/DCM_mono_poly/data/UKBB_DCM_mutations_and_score_v3.tsv")
ukbb <- as.data.frame(ukbb)[, c(
  "IID", "DCM", "Heart_Failure", "Hypertrophic_cardiomyopathy", "age", "sex",
  "exclude_sampleQC", "exclude_related_2nd_deg", "in_white_British_ancestry_subset",
  "LDL_direct",
  "f.21001.0.0",  # BMI
  "f.4080.0.0", "f.4079.0.0",  # Blood pressure
  "Coronary_Artery_Disease_INTERMEDIATE", "Diabetes_Type_2",
  "Atrial_fibrillation_or_flutter", "Chronic_kidney_disease", "Ischemic_stroke"
)]

# =======================
# 3. Map app7089 IDs to app17488 IDs
# =======================
linker <- fread("/mnt/project/sjj/short_projects/DCM_mono_poly/data/ukb_app17488_app7089_link.csv")
ukbb <- merge(ukbb, linker, by.x = "IID", by.y = "app7089", all.x = TRUE)

# =======================
# 4. Assign Exome Sequencing Batches
# =======================
ukbb$batch <- "450k"
u200k <- fread("/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_200k_tranche.fam", header = FALSE)
u50k  <- fread("/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_50k_tranche.fam",  header = FALSE)

ukbb$batch[ukbb$app17488 %in% u200k$V1] <- "200k"
ukbb$batch[ukbb$IID %in% u50k$V1] <- "50k"

# =======================
# 5. Merge Race/Ethnicity Information
# =======================
race <- fread("/mnt/project/kbiddinger/data/UKBBRace.txt")
ukbb <- merge(ukbb, race[, c("eid", setdiff(colnames(race), colnames(ukbb)))], 
              by.x = "IID", by.y = "eid", all.x = TRUE)

# =======================
# 6. Clean and Rename Variables
# =======================
colnames(ukbb)[colnames(ukbb) == "LDL_direct"] <- "LDL_C"
colnames(ukbb)[colnames(ukbb) == "f.21001.0.0"] <- "BMI"
colnames(ukbb)[colnames(ukbb) == "f.4080.0.0"] <- "SBP"
colnames(ukbb)[colnames(ukbb) == "f.4079.0.0"] <- "DBP"
colnames(ukbb)[colnames(ukbb) == "Coronary_Artery_Disease_INTERMEDIATE"] <- "CAD"
colnames(ukbb)[colnames(ukbb) == "Diabetes_Type_2"] <- "T2D"
colnames(ukbb)[colnames(ukbb) == "Atrial_fibrillation_or_flutter"] <- "AF"
colnames(ukbb)[colnames(ukbb) == "Chronic_kidney_disease"] <- "CKD"
colnames(ukbb)[colnames(ukbb) == "Ischemic_stroke"] <- "IS"

# =======================
# 7. Apply Population-Level QC
# =======================
ukbb <- subset(ukbb, exclude_sampleQC == 0)
ukbb <- subset(ukbb, is.na(exclude_related_2nd_deg))

# =======================
# 8. Format IDs for REGENIE Compatibility
# =======================
ukbb$FID <- ukbb$app17488
ukbb$IID <- ukbb$app17488

# =======================
# 9. Load and Standardize PCA Components
# =======================
PCs_to_remove <- paste0("PC", 1:10)
ukbb <- ukbb[, !(names(ukbb) %in% PCs_to_remove)]

pca <- fread("/mnt/project/exome_450k_plink/PCA/PCs_round1.txt")
pca$sample_id <- gsub("_.*", "", pca$sample_id)
ukbb <- merge(ukbb, pca, by.x = "IID", by.y = "sample_id")

PCs <- paste0("PC", 1:20)
ukbb[PCs] <- scale(ukbb[PCs])

# =======================
# 10. Add Derived Covariates
# =======================
ukbb$age_squared <- ukbb$age^2

# =======================
# 11. Define Traits of Interest
# =======================
qt <- c("LDL_C", "BMI", "SBP", "DBP")
bt <- c("CAD", "T2D", "AF", "CKD", "IS")

# =======================
# 12. Remove Outliers from Quantitative Traits
# =======================
for (t in qt) {
  vec <- ukbb[[t]]
  iqr <- IQR(vec, na.rm = TRUE)
  upper <- quantile(vec, 0.75, na.rm = TRUE) + 1.5 * iqr
  lower <- quantile(vec, 0.25, na.rm = TRUE) - 1.5 * iqr
  ukbb[[t]][vec > upper | vec < lower] <- NA
}

# =======================
# 13. Save and Upload Histograms
# =======================
for (t in qt) {
  png(paste0(t, ".png"))
  hist(
    ukbb[[t]],
    main = "",
    xlab = "",
    cex.lab = 1.5,
    cex.axis = 1.6
  )
  title(xlab = t, font.lab = 2, cex.lab = 2)
  dev.off()
  
  system(paste0("dx rm exome-seq:/kbiddinger/projects/mashHF/data/pheno/phewas/", t, ".png"))
  system(paste0("dx upload ", t, ".png --destination exome-seq:kbiddinger/projects/mashHF/data/pheno/phewas/"))
}

# =======================
# 14. Export Final Curated Dataset
# =======================
file_name <- "mashHF_pheno_phewas.tsv"
write.table(ukbb[, c("FID", "IID", qt, bt, "age", "age_squared", "sex", PCs, "batch")],
            file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)

system(paste0("dx rm exome-seq:/kbiddinger/projects/mashHF/data/pheno/phewas/", file_name))
system(paste0("dx upload ", file_name, " --destination exome-seq:kbiddinger/projects/mashHF/data/pheno/phewas/"))
