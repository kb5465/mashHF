# Load required libraries
install.packages('R.utils')
library('data.table')
install.packages('bit64')
library(dplyr)

# Read in bulk data file
ukbb <- fread("/mnt/project/sjj/short_projects/DCM_mono_poly/data/UKBB_DCM_mutations_and_score_v3.tsv")
ukbb <- as.data.frame(ukbb)[,c("IID", "DCM", "Heart_Failure", "Hypertrophic_cardiomyopathy", "age", "sex",
                               "exclude_sampleQC", "exclude_related_2nd_deg", "in_white_British_ancestry_subset",
                               "LDL_direct",
                               "f.21001.0.0", # BMI
                               "f.4080.0.0", "f.4079.0.0", # Blood pressure
                               "Coronary_Artery_Disease_INTERMEDIATE", "Diabetes_Type_2", 
                               "Atrial_fibrillation_or_flutter", "Chronic_kidney_disease",
                               "Ischemic_stroke"
                               )]

# Load linker file
linker <- fread('/mnt/project/sjj/short_projects/DCM_mono_poly/data/ukb_app17488_app7089_link.csv')
ukbb <- merge(ukbb, linker, by.x = "IID", by.y = "app7089", all.x = TRUE)

# Get exome sequencing batch
ukbb$batch <- "450k"
u200k <- fread('/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_200k_tranche.fam',
               stringsAsFactors = F, data.table = F, header = F)
ukbb$batch[ukbb$app17488 %in% u200k$V1] <- "200k"
u50k <- fread('/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_50k_tranche.fam',
              stringsAsFactors = F, data.table = F, header = F)
ukbb$batch[ukbb$IID %in% u50k$V1] <- "50k"

# Merge race variable
race <- as.data.frame(fread("/mnt/project/kbiddinger/data/UKBBRace.txt"))
ukbb <- merge(ukbb,race[,c("eid",setdiff(colnames(race),colnames(ukbb)))], by.x = 'IID', by.y = 'eid', all.x = TRUE)

# Rename variables
colnames(ukbb)[which(colnames(ukbb) == "LDL_direct")] <- "LDL_C"
colnames(ukbb)[which(colnames(ukbb) == "f.21001.0.0")] <- "BMI"
colnames(ukbb)[which(colnames(ukbb) == "f.4080.0.0")] <- "SBP"
colnames(ukbb)[which(colnames(ukbb) == "f.4079.0.0")] <- "DBP"
colnames(ukbb)[which(colnames(ukbb) == "Coronary_Artery_Disease_INTERMEDIATE")] <- "CAD"
colnames(ukbb)[which(colnames(ukbb) == "Diabetes_Type_2")] <- "T2D"
colnames(ukbb)[which(colnames(ukbb) == "Atrial_fibrillation_or_flutter")] <- "AF"
colnames(ukbb)[which(colnames(ukbb) == "Chronic_kidney_disease")] <- "CKD"
colnames(ukbb)[which(colnames(ukbb) == "Ischemic_stroke")] <- "IS"

# Population QC using centralized UKBB metrics
ukbb <- subset(ukbb, ukbb$exclude_sampleQC==0)
ukbb <- subset(ukbb, is.na(ukbb$exclude_related_2nd_deg))

# Format ID variable for REGENIE software
ukbb$FID <- ukbb$app17488
ukbb$IID <- ukbb$app17488

# Remove old PCs and load new PCs
PCs <- paste0("PC", seq(1,10))
ukbb <- ukbb[, !(names(ukbb) %in% PCs)]
pca <- fread('/mnt/project/exome_450k_plink/PCA/PCs_round1.txt', stringsAsFactors = F, data.table = F)
pca$sample_id <- gsub("_.*", "", pca$sample_id)
ukbb <- merge(ukbb, pca, by.x = 'IID', by.y = 'sample_id')
PCs <- paste0("PC", seq(1,20))
ukbb[PCs] <- scale(ukbb[PCs])

# Construct age squared for covariates
ukbb$age_squared <- ukbb$age ^ 2

# Define outcomes of interest
qt <- c("LDL_C", "BMI", "SBP", "DBP")
bt <- c("CAD", "T2D", "AF", "CKD", "IS")

# Remove outliers
for (t in qt) {
  # Pull column
  vec <- ukbb[[t]]
  
  # Calculate IQR
  iqr <- IQR(vec, na.rm = T)
  
  # Calculate upper and lower bounds of outliers
  upper <- quantile(vec, probs=.75, na.rm = T) + 1.5 * iqr
  lower <- quantile(vec, probs=.25, na.rm = T) - 1.5 * iqr
  
  # Set outliers to NA
  ukbb[[t]][vec > upper | vec < lower] <- NA
}

# Inspect data distributions
for (t in qt) {
  output_name <- paste0(t, ".png")
  
  # Save the plot as a PNG
  png(output_name)
  
  # Customize the plot
  hist(
    ukbb[[t]],
    main = "",
    xlab = "",
    cex.lab = 1.5,
    cex.main = 1.8,
    cex.axis = 1.6
  )
  
  # Make axis titles bold
  title(xlab = t, font.lab = 2, cex.lab = 2)  # Bold x-axis label
  
  dev.off()
  
  system(paste0("dx rm exome-seq:/kbiddinger/projects/mashHF/data/pheno/phewas/", output_name))
  system(paste0("dx upload ", output_name, " --destination exome-seq:kbiddinger/projects/mashHF/data/pheno/phewas/"))
}

# Write the final curated data to a tab-separated file
file_name <- "mashHF_pheno_phewas.tsv"
write.table(ukbb[, c("FID", "IID", qt, bt,
                     "age", "age_squared", "sex", PCs, "batch")], file_name, 
            sep = "\t", col.names = T, row.names = F, quote = F)

# Upload the curated data file to a destination
system(paste0("dx rm exome-seq:/kbiddinger/projects/mashHF/data/pheno/phewas/", file_name))
system(paste0("dx upload ", file_name, " --destination exome-seq:kbiddinger/projects/mashHF/data/pheno/phewas/"))
