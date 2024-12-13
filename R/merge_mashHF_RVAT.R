library(data.table)
install.packages("mashr")
library(mashr)

# File directory
folder_path <- "/mnt/project/kbiddinger/projects/mashHF/results/prior/"

# List files in folder
file_list <- list.files(folder_path, full.names = TRUE)
file_list <- subset(file_list, !grepl("LVConc", file_list))

# Empty list for tables
data_list <- list()

# Loop files and read into tables
for (file in file_list) {
  # Read the file
  dt <- fread(file)
  
  # Add the table to the list
  data_list[[file]] <- dt
}

# Merge all tables by ID column
final_merged <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_list)

# Print TTN gene row
subset(final_merged, final_merged$ID == "ENSG00000155657")

# Write the final curated data to a tab-separated file
file_name <- "mashHF_canonical_LOF_0.01_merged.tsv"
write.table(final_merged, file_name, sep = "\t", col.names = T, row.names = F, quote = F)

# Upload the curated data file to a destination
system(paste0("dx rm exome-seq:/kbiddinger/projects/mashHF/data/", file_name))
system(paste0("dx upload ", file_name, " --destination exome-seq:kbiddinger/projects/mashHF/data/"))
