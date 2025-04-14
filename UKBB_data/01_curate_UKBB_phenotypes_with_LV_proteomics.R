###############################################
# mashHF: UKBB Phenotype Preparation Script
# Author: Kiran Biddinger
# Description: Load, merge, and clean phenotype, imaging,
# and proteomics data from the UK Biobank for mashHF analysis.
###############################################

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# =======================
# 1. Load Core Phenotype Data
# =======================

# Read main phenotype and genetics data
ukbb <- fread("/mnt/project/sjj/short_projects/DCM_mono_poly/data/UKBB_DCM_mutations_and_score_v3.tsv")
ukbb <- as.data.frame(ukbb)[, c(
  "IID", "DCM", "Heart_Failure", "Hypertrophic_cardiomyopathy", "age", "sex",
  "exclude_sampleQC", "exclude_related_2nd_deg", "in_white_British_ancestry_subset",
  "f.21001.0.0", "f.20116.0.0", "LDL_direct",
  "f.4080.0.0", "f.4079.0.0",  # Blood pressure
  "Coronary_Artery_Disease_INTERMEDIATE", "Diabetes_Type_2",
  "Atrial_fibrillation_or_flutter", "Chronic_kidney_disease", "Ischemic_stroke"
)]

# =======================
# 2. Merge Identifiers & Sequencing Batches
# =======================

# Map app7089 to app17488 IDs
linker <- fread("/mnt/project/sjj/short_projects/DCM_mono_poly/data/ukb_app17488_app7089_link.csv")
ukbb <- merge(ukbb, linker, by.x = "IID", by.y = "app7089", all.x = TRUE)

# Assign sequencing batch
ukbb$batch <- "450k"
u200k <- fread('/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_200k_tranche.fam')
u50k  <- fread('/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_50k_tranche.fam')
ukbb$batch[ukbb$app17488 %in% u200k$V1] <- "200k"
ukbb$batch[ukbb$IID %in% u50k$V1] <- "50k"

# =======================
# 3. Merge Proteomics (NTproBNP, Troponin)
# =======================

load_olink <- function(filepath, prefix) {
  df <- fread(filepath)
  df <- as.data.frame(df)
  colnames(df) <- c("eid", paste0("ins_index_", prefix), paste0("protein_id_", prefix),
                    prefix, paste0("Batch_", prefix), paste0("plate_", prefix))
  df[!duplicated(df$eid), ]
}

NTproBNP <- load_olink("/mnt/project/kbiddinger/data/olink/ukbb_proteomics_NTproBNP.tsv", "NTproBNP")
Troponin <- load_olink("/mnt/project/kbiddinger/data/olink/ukbb_proteomics_Troponin.tsv", "Troponin")

ukbb <- merge(ukbb, NTproBNP, by.x = "app17488", by.y = "eid", all.x = TRUE)
ukbb <- merge(ukbb, Troponin, by.x = "app17488", by.y = "eid", all.x = TRUE)

# =======================
# 4. Merge Cardiac MRI Traits and Derived Metrics
# =======================

cMRI <- fread("/mnt/project/kbiddinger/projects/mashHF/data/CMRI.csv")
bsa  <- fread("/mnt/project/kbiddinger/projects/mashHF/data/UKBB_BSA.tsv")
lvm  <- fread("/mnt/project/kbiddinger/projects/mashHF/data/seg_lvm_age.csv")

ukbb <- merge(ukbb, cMRI[, c("app17488", "LVESV", "LVEDV", "LVEF",
                             "LV_mean_myocardial_wall_thickness_global", 
                             "LV_circumferential_strain_global", "LV_longitudinal_strain_global", 
                             "LV_radial_strain_global")], 
              by = "app17488", all.x = TRUE)

ukbb <- merge(ukbb, lvm[, c("sample_id", "lvm_seg_adjusted", "lvmi_seg_adjusted")], 
              by.x = "IID", by.y = "sample_id", all.x = TRUE)

ukbb <- merge(ukbb, bsa, by.x = "IID", by.y = "f.eid", all.x = TRUE)

# =======================
# 5. Add Race and Derived Variables
# =======================

race <- fread("/mnt/project/kbiddinger/data/UKBBRace.txt")
ukbb <- merge(ukbb, race[, c("eid", setdiff(colnames(race), colnames(ukbb)))], 
              by.x = "IID", by.y = "eid", all.x = TRUE)

# Define derived cardiac traits
ukbb$LVConc <- ukbb$LVEDV / ukbb$lvm_seg_adjusted
ukbb$LVESVi <- ukbb$LVESV / ukbb$bsa
ukbb$LVEDVi <- ukbb$LVEDV / ukbb$bsa

# Rename columns
colnames(ukbb)[colnames(ukbb) == "lvmi_seg_adjusted"] <- "LVMi"
colnames(ukbb)[colnames(ukbb) == "LV_mean_myocardial_wall_thickness_global"] <- "meanWT"
colnames(ukbb)[colnames(ukbb) == "LV_circumferential_strain_global"] <- "strain_circ"
colnames(ukbb)[colnames(ukbb) == "LV_longitudinal_strain_global"] <- "strain_long"
colnames(ukbb)[colnames(ukbb) == "LV_radial_strain_global"] <- "strain_rad"
colnames(ukbb)[colnames(ukbb) == "Hypertrophic_cardiomyopathy"] <- "HCM"
colnames(ukbb)[colnames(ukbb) == "Heart_Failure"] <- "HF"

# =======================
# 6. Quality Control and Final Covariates
# =======================

# Filter on centralized QC
ukbb <- subset(ukbb, exclude_sampleQC == 0)
ukbb <- subset(ukbb, is.na(exclude_related_2nd_deg))

# Format IDs for REGENIE
ukbb$FID <- ukbb$app17488
ukbb$IID <- ukbb$app17488

# Merge PCA components and standardize
PCs_to_remove <- paste0("PC", 1:10)
ukbb <- ukbb[, !(names(ukbb) %in% PCs_to_remove)]
pca <- fread('/mnt/project/exome_450k_plink/PCA/PCs_round1.txt')
pca$sample_id <- gsub("_.*", "", pca$sample_id)
ukbb <- merge(ukbb, pca, by.x = "IID", by.y = "sample_id")
ukbb[paste0("PC", 1:20)] <- scale(ukbb[paste0("PC", 1:20)])

# Add age squared term
ukbb$age_squared <- ukbb$age^2

# =======================
# 7. Define and Clean Traits
# =======================

qt <- c("LVMi", "LVESVi", "LVEDVi", "LVEF", "meanWT",
        "strain_rad", "strain_long", "strain_circ", "Troponin", "NTproBNP")
bt <- c("HF", "DCM", "HCM")

# Remove outliers using IQR
for (t in qt) {
  vec <- ukbb[[t]]
  iqr <- IQR(vec, na.rm = TRUE)
  upper <- quantile(vec, 0.75, na.rm = TRUE) + 1.5 * iqr
  lower <- quantile(vec, 0.25, na.rm = TRUE) - 1.5 * iqr
  ukbb[[t]][vec > upper | vec < lower] <- NA
}

# =======================
# 8. Save and Upload Final Curated Dataset
# =======================

file_name <- "mashHF_pheno.tsv"
write.table(ukbb[, c("FID", "IID", qt, bt, "age", "age_squared", "sex", paste0("PC", 1:20), "batch")], 
            file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
system(paste0("dx rm exome-seq:/kbiddinger/projects/mashHF/data/pheno/", file_name))
system(paste0("dx upload ", file_name, " --destination exome-seq:kbiddinger/projects/mashHF/data/pheno/"))

# =======================
# 9. Density Plots: 3x3 + 1 Layout
# =======================

# Generate density plots
plot_list <- lapply(qt, function(t) {
  vec <- ukbb[[t]]
  vec <- vec[!is.na(vec)]
  if (length(vec) < 10) return(NULL)
  ggplot(data.frame(val = vec), aes(x = val)) +
    geom_density(fill = "steelblue", alpha = 0.3, color = "black") +
    labs(title = t, x = t, y = "Density") +
    theme_minimal(base_size = 16) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10))
})

# Filter valid plots
plots <- plot_list[!sapply(plot_list, is.null)]

# Create 3x3 grid and middle-bottom layout
top_grid <- ggarrange(plotlist = plots[1:9], ncol = 3, nrow = 3)
bottom_row <- ggarrange(ggplot() + theme_void(), plots[[10]], ggplot() + theme_void(), ncol = 3)
final_plot <- ggarrange(top_grid, bottom_row, nrow = 2, heights = c(3, 1))

# Save and upload
ggsave("combined_density_plots.png", final_plot, width = 15, height = 18, dpi = 300)
system("dx rm exome-seq:/kbiddinger/projects/mashHF/data/pheno/combined_density_plots.png")
system("dx upload combined_density_plots.png --destination exome-seq:kbiddinger/projects/mashHF/data/pheno/")
