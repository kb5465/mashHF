#!/usr/bin/env Rscript

###############################################
# mashHF: Merge Canonical LOF 0.01 Gene-Level Results (Firth)
# Author: Kiran Biddinger
# Description:
#   - Reads all gene-level result files from `st/` folder
#   - Merges them by gene ID (ENSG)
#   - Excludes LVConc to avoid issues
#   - Outputs a single TSV file and uploads to DNAnexus
###############################################

# ----------------------------
# 1. Load Required Libraries
# ----------------------------
library(data.table)

# ----------------------------
# 2. Define Input Folder
# ----------------------------
folder_path <- "/mnt/project/kbiddinger/projects/mashHF/results/st/"

# Get list of all result files, excluding any with "LVConc"
file_list <- list.files(folder_path, full.names = TRUE)
file_list <- subset(file_list, !grepl("LVConc", file_list))

# ----------------------------
# 3. Read All Files into a List
# ----------------------------
data_list <- list()
for (file in file_list) {
  dt <- fread(file)
  data_list[[file]] <- dt
}

# ----------------------------
# 4. Merge All Tables by Gene ID
# ----------------------------
final_merged <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_list)

# ----------------------------
# 5. Print TTN Gene Row for Inspection
# ----------------------------
ttn_row <- subset(final_merged, ID == "ENSG00000155657")
print(ttn_row)

# ----------------------------
# 6. Write Final Table to TSV File
# ----------------------------
file_name <- "mashHF_canonical_LOF_0.01_merged_firth.tsv"
write.table(final_merged, file_name, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# ----------------------------
# 7. Upload Output File to DNAnexus
# ----------------------------
system(paste0("dx rm exome-seq:/kbiddinger/projects/mashHF/data/", file_name))
system(paste0("dx upload ", file_name, " --destination exome-seq:kbiddinger/projects/mashHF/data/"))
