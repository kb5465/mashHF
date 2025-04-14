###############################################
# mashHF: Table 1 Summary Statistics
# Author: Kiran Biddinger
# Description: Generate descriptive statistics
# for the full cohort and key analysis subsets
# using the tableone package.
###############################################

# =======================
# 1. Load Required Library
# =======================
# Only needed once per environment; skip if already installed
# install.packages('tableone') 
library(tableone)

# =======================
# 2. Define Subsets
# =======================
# Subset for individuals with LV ejection fraction data
ukbb$lv_subset <- !is.na(ukbb$LVEF)

# Subset for individuals with machine learning-derived LV mass
ukbb$lv_ml_subset <- !is.na(ukbb$LVMi)

# Subset for individuals with proteomics data
ukbb$prot_subset <- !is.na(ukbb$NTproBNP)

# =======================
# 3. Define Variables for Summary
# =======================
# Variables to summarize
myvars <- c("sex", "age", "Race", 
            "f.21001.0.0.x",  # BMI
            "f.20116.0.0.x",  # Smoking status
            "f.4080.0.0.x", "f.4079.0.0.x",  # SBP, DBP
            "LDL_direct", "LVEF",
            "Coronary_Artery_Disease_INTERMEDIATE", "Diabetes_Type_2",
            "Atrial_fibrillation_or_flutter", "Chronic_kidney_disease",
            "Ischemic_stroke", "HF", "DCM", "HCM")

# Categorical variables for table stratification
catvars <- c("sex", "Race", "f.20116.0.0.x",
             "Coronary_Artery_Disease_INTERMEDIATE", "Diabetes_Type_2",
             "Atrial_fibrillation_or_flutter", "Chronic_kidney_disease",
             "Ischemic_stroke", "HF", "DCM", "HCM")

# =======================
# 4. Create Subset Grouping Variable
# =======================
ukbb$subset_group <- "All Participants"
ukbb$subset_group[ukbb$lv_subset]     <- "LV Subset"
ukbb$subset_group[ukbb$lv_ml_subset]  <- "LV ML Subset"
ukbb$subset_group[ukbb$prot_subset]   <- "Proteomics Subset"

# Set factor levels in the desired order
ukbb$subset_group <- factor(
  ukbb$subset_group,
  levels = c("All Participants", "LV Subset", "LV ML Subset", "Proteomics Subset")
)

# =======================
# 5. Create Table 1 for Each Subset
# =======================
options(digits = 1)  # Set number of digits for output

tab_all   <- print(CreateTableOne(vars = myvars, data = ukbb, factorVars = catvars), printToggle = FALSE)
tab_lv    <- print(CreateTableOne(vars = myvars, data = ukbb[ukbb$lv_subset, ], factorVars = catvars), printToggle = FALSE)
tab_lvml  <- print(CreateTableOne(vars = myvars, data = ukbb[ukbb$lv_ml_subset, ], factorVars = catvars), printToggle = FALSE)
tab_prot  <- print(CreateTableOne(vars = myvars, data = ukbb[ukbb$prot_subset, ], factorVars = catvars), printToggle = FALSE)

# =======================
# 6. Combine and Export Table
# =======================
# Combine tables into one CSV
combined_tab <- cbind(tab_all, tab_lv, tab_lvml, tab_prot)
write.csv(combined_tab, "UKBB.MashHF.Table1.FULL.csv")

# Upload to project folder on DNAnexus
system("dx upload UKBB.MashHF.Table1.FULL.csv --destination exome-seq:kbiddinger/projects/mashHF/results/")
