install.packages('tableone')
library(tableone)

# Create subsets
ukbb$lv_subset <- !is.na(ukbb$LVEF)
ukbb$lv_ml_subset <- !is.na(ukbb$LVMi)
ukbb$prot_subset <- !is.na(ukbb$NTproBNP)

# Define variables
myvars <- c("sex","age","Race", "f.21001.0.0.x", "f.20116.0.0.x", "f.4080.0.0.x", "f.4079.0.0.x", "LDL_direct",
            "LVEF",
            "Coronary_Artery_Disease_INTERMEDIATE", "Diabetes_Type_2", 
            "Atrial_fibrillation_or_flutter", "Chronic_kidney_disease",
            "Ischemic_stroke",
            "HF", "DCM", "HCM")
catvars <- c("sex","Race", "f.20116.0.0.x",  
             "Coronary_Artery_Disease_INTERMEDIATE", "Diabetes_Type_2", 
             "Atrial_fibrillation_or_flutter", "Chronic_kidney_disease",
             "Ischemic_stroke",
             "HF", "DCM", "HCM")

# Create a new grouping variable
ukbb$subset_group <- "All Participants"
ukbb$subset_group[ukbb$lv_subset] <- "LV Subset"
ukbb$subset_group[ukbb$lv_ml_subset] <- "LV ML Subset"
ukbb$subset_group[ukbb$prot_subset] <- "Proteomics Subset"

# Make sure the group is a factor with the desired order
ukbb$subset_group <- factor(ukbb$subset_group, levels = c("All Participants", "LV Subset", "LV ML Subset", "Proteomics Subset"))

# Create Table 1
options(digits=1)
tab_all <- print(CreateTableOne(vars = myvars, data = ukbb, factorVars = catvars), printToggle = FALSE)
tab_lv <- print(CreateTableOne(vars = myvars, data = ukbb[ukbb$lv_subset,], factorVars = catvars), printToggle = FALSE)
tab_lvml <- print(CreateTableOne(vars = myvars, data = ukbb[ukbb$lv_ml_subset,], factorVars = catvars), printToggle = FALSE)
tab_prot <- print(CreateTableOne(vars = myvars, data = ukbb[ukbb$prot_subset,], factorVars = catvars), printToggle = FALSE)

# Combine them
combined_tab <- cbind(tab_all, tab_lv, tab_lvml, tab_prot)
write.csv(combined_tab, "UKBB.MashHF.Table1.FULL.csv")
system(paste0("dx upload ", "UKBB.MashHF.Table1.FULL.csv", " --destination exome-seq:kbiddinger/projects/mashHF/results/"))

