# Load required libraries
install.packages('R.utils')
library('data.table')
install.packages('bit64')
library(dplyr)

# Read in bulk data file
ukbb <- fread("/mnt/project/sjj/short_projects/DCM_mono_poly/data/UKBB_DCM_mutations_and_score_v3.tsv")
ukbb <- as.data.frame(ukbb)[,c("IID", "DCM", "Heart_Failure", "Hypertrophic_cardiomyopathy", "age", "sex",
                               "exclude_sampleQC", "exclude_related_2nd_deg")]

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

# Load olink proteomics data
NTproBNP <- fread("/mnt/project/kbiddinger/data/olink/ukbb_proteomics_NTproBNP.tsv")
NTproBNP <- as.data.frame(NTproBNP)
colnames(NTproBNP) <- c("eid", "ins_index_NTproBNP","protein_id_NTproBNP", "NTproBNP", "Batch_NTproBNP", "plate_NTproBNP")
NTproBNP <- NTproBNP[!duplicated(NTproBNP$eid),]
Troponin <- fread("/mnt/project/kbiddinger/data/olink/ukbb_proteomics_Troponin.tsv")
Troponin <- as.data.frame(Troponin)
colnames(Troponin) <- c("eid", "ins_index_Troponin","protein_id_Troponin", "Troponin", "Batch_Troponin", "plate_Troponin")
Troponin <- Troponin[!duplicated(Troponin$eid),]
ukbb <- merge(ukbb, NTproBNP, by.x = "app17488", by.y = "eid", all.x=T)
ukbb <- merge(ukbb, Troponin, by.x = "app17488", by.y = "eid", all.x=T)

# Load LV trait data
cMRI <- fread("/mnt/project/kbiddinger/projects/mashHF/data/CMRI.csv")
bsa <- fread("/mnt/project/kbiddinger/projects/mashHF/data/UKBB_BSA.tsv")
lvm <- fread("/mnt/project/kbiddinger/projects/mashHF/data/seg_lvm_age.csv")

# Extract cMRI parameters
ukbb <- merge(ukbb, as.data.frame(cMRI)[,c("app17488", "LVESV", "LVEDV", "LVEF",
                                           "LV_mean_myocardial_wall_thickness_global", 
                                           "LV_circumferential_strain_global", "LV_longitudinal_strain_global", "LV_radial_strain_global")], 
              by.x = "app17488", by.y = "app17488", all.x=T)
ukbb <- merge(ukbb, as.data.frame(lvm)[,c("sample_id", "lvm_seg_adjusted", "lvmi_seg_adjusted")], by.x = "IID", by.y = "sample_id", all.x=T)
ukbb <- merge(ukbb, bsa, by.x = "IID", by.y = "f.eid", all.x=T)

# Define LV concentricity
ukbb$LVConc <- ukbb$LVEDV / ukbb$lvm_seg_adjusted

# Index LV volumes to body surface area
ukbb$LVESVi <- ukbb$LVESV / ukbb$bsa
ukbb$LVEDVi <- ukbb$LVEDV / ukbb$bsa

# Rename variables
colnames(ukbb)[which(colnames(ukbb) == "lvmi_seg_adjusted")] <- "LVMi"
colnames(ukbb)[which(colnames(ukbb) == "LV_mean_myocardial_wall_thickness_global")] <- "meanWT"
colnames(ukbb)[which(colnames(ukbb) == "LV_circumferential_strain_global")] <- "strain_circ"
colnames(ukbb)[which(colnames(ukbb) == "LV_longitudinal_strain_global")] <- "strain_long"
colnames(ukbb)[which(colnames(ukbb) == "LV_radial_strain_global")] <- "strain_rad"
colnames(ukbb)[which(colnames(ukbb) == "Hypertrophic_cardiomyopathy")] <- "HCM"
colnames(ukbb)[which(colnames(ukbb) == "Heart_Failure")] <- "HF"

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
qt <- c("LVMi", "LVESVi", "LVEDVi", "LVEF", "LVConc",
             "meanWT", "strain_rad", "strain_long", "strain_circ", "Troponin", "NTproBNP")
bt <- c("HF", "DCM", "HCM")

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
  hist(ukbb[[t]], main = t)
}

# Write the final curated data to a tab-separated file
file_name <- "mashHF_pheno.tsv"
write.table(ukbb[, c("FID", "IID", qt, bt,
                     "age", "age_squared", "sex", PCs, "batch")], file_name, 
            sep = "\t", col.names = T, row.names = F, quote = F)

# Upload the curated data file to a destination
system(paste0("dx rm exome-seq:/kbiddinger/projects/mashHF/data/pheno/", file_name))
system(paste0("dx upload ", file_name, " --destination exome-seq:kbiddinger/projects/mashHF/data/pheno/"))
