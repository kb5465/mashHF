# ========================================
# Title: Canonical LOF Gene-Trait Heatmap
# Description: Loads PheWAS results, computes Z-scores, and visualizes associations
# Author: Kiran Biddinger
# ========================================

# ---- Load Libraries ----
library(tidyverse)
library(data.table)
library(biomaRt)

# ---- File Settings ----
# Path to directory with TSVs
data_dir <- "/Users/kbidding/Documents/mashHF/data/pheWAS"

# Get list of all canonical_LOF result files
files <- list.files(data_dir, pattern = "_canonical_LOF_0.01_merged.tsv$", full.names = TRUE)

# Extract phenotype names (PHENO)
phenos <- gsub("_canonical_LOF_0.01_merged.tsv", "", basename(files))
names(files) <- phenos

# ---- Load All Files and Compute Z-scores ----
# Function to read file and compute z-scores
read_and_zscore <- function(file, pheno) {
  df <- as.data.frame(fread(file))
  
  # Compute Z-scores as beta / se
  df[[pheno]] <- df[[paste0("BETA_", pheno)]] / df[[paste0("SE_", pheno)]]
  
  # Keep only gene ID and Z-score
  df_subset <- df[, c("ID", pheno)]
  return(df_subset)
}

# Apply the function to all files
zscore_list <- map2(files, phenos, read_and_zscore)

# Merge all data frames by gene ID
zscore_data <- reduce(zscore_list, full_join, by = "ID")

# ---- Annotate Genes with Gene Names ----
# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch gene annotations
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = zscore_data$ID,
  mart = ensembl
)

# Rename for merging
colnames(gene_annotations) <- c("ID", "Gene_Name")

# Merge to get gene names
zscore_data <- left_join(zscore_data, gene_annotations, by = "ID")

# Reorder columns: Gene_Name first
zscore_data <- zscore_data %>%
  select(Gene_Name, all_of(phenos))

genes_heatmap_fn <- function(data, plotname) {
  # Convert to long format
  heatmap_data <- data %>%
    pivot_longer(cols = -Gene_Name, names_to = "Trait", values_to = "Z")
  
  # Compute Bonferroni-corrected Z threshold
  alpha <- 0.05 / 99
  z_thresh <- qnorm(1 - alpha / 2)
  
  # Set max Z for color scaling (symmetric)
  max_z <- ceiling(max(abs(heatmap_data$Z), na.rm = TRUE))
  
  # Create breaks for color scale (evenly spaced and including ±z_thresh)
  breaks <- sort(unique(c(-max_z, -z_thresh, -round(max_z/2, 1), 0, round(max_z/2, 1), z_thresh, max_z)))
  labels <- format(round(breaks, 1), nsmall = 0, trim = TRUE)
  
  # Factor levels for plotting
  heatmap_data$Gene_Name <- factor(heatmap_data$Gene_Name, levels = rev(sort(unique(heatmap_data$Gene_Name))))
  heatmap_data$Trait <- factor(heatmap_data$Trait, levels = c("SBP", "DBP", "LDL_C", "BMI", "AF", "CKD", "CAD", "IS", "T2D"))
  
  # Plot
  p <- ggplot(heatmap_data, aes(x = Trait, y = Gene_Name, fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "red", mid = "white", high = "blue", midpoint = 0,
      limits = c(-max_z, max_z),
      breaks = breaks,
      labels = labels,
      name = "Z-score"
    ) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"))
  
  ggsave(plotname, p, width = 6.56, height = 12)
  print(p)
}

# ---- Run and Save Heatmap ----
genes_heatmap_fn(
  zscore_data,
  plotname = "/Users/kbidding/Documents/mashHF/Figures/gene_trait_zscore_heatmap.pdf"
)
