#!/usr/bin/env Rscript

################################################################################
# Title: Merge Single-Trait RVAT Results
# Description:
#   Merges REGENIE Step 2 results across multiple phenotypes for canonical LOF
#   genes (MAF < 1%). Produces a single merged table by Ensembl ID and uploads
#   the file to a DNAnexus project folder.
#
# Author: Kiran Biddinger
# Date: April 2025
################################################################################

# -----------------------------
# Load required libraries
# -----------------------------
library(data.table)

# -----------------------------
# Set file directory
# -----------------------------
folder_path <- "/mnt/project/kbiddinger/projects/mashHF/results/prior/"

# -----------------------------
# List and filter result files
# -----------------------------
file_list <- list.files(folder_path, full.names = TRUE)
file_list <- file_list[!grepl("LVConc", file_list)]  # Exclude known problematic file

# -----------------------------
# Read and store all files
# -----------------------------
data_list <- lapply(file_list, fread)

# -----------------------------
# Merge all tables by 'ID'
# -----------------------------
final_merged <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_list)

# -----------------------------
# (Optional) Inspect specific genes
# -----------------------------
print(subset(final_merged, ID == "ENSG00000155657"))  # TTN

# -----------------------------
# Save merged table
# -----------------------------
output_file <- "mashHF_canonical_LOF_0.01_merged_firth.tsv"
write.table(final_merged, file = output_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# -----------------------------
# Upload to DNAnexus
# -----------------------------
# Remove any existing file before upload (safe cleanup)
system(paste0("dx rm -f exome-seq:/kbiddinger/projects/mashHF/data/", output_file))
system(paste0("dx upload ", output_file, " --destination exome-seq:/kbiddinger/projects/mashHF/data/"))
