#!/usr/bin/env Rscript

###############################################################
# mashHF: Heatmap of Gene-Trait Associations (Z-scores)
# Author: Kiran Biddinger
# Description:
# - Identifies genes significantly associated (Bonferroni) with ≥1 trait
# - Plots a heatmap of Z-scores across traits
# - Optional 2-panel display separating disease vs. quantitative traits
###############################################################

# ----------------------------
# Load Required Libraries
# ----------------------------
library(ggplot2)
library(dplyr)
library(data.table)
library(tibble)
library(tidyr)
library(RColorBrewer)

# ----------------------------
# Helper: Identify Significant Genes (p < 0.05/m for any trait)
# ----------------------------
get_all_genes <- function(data) {
  results <- c()
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    z <- data[[paste0("Z_", pheno)]]
    pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
    bonf <- 0.05 / sum(!is.na(pval))
    genes <- data$Gene_Name[pval < bonf]
    results <- c(results, genes)
  }
  results <- unique(results[results != ""])
  return(results)
}

# ----------------------------
# Helper: Extract Z-scores for Significant Genes
# ----------------------------
get_gene_associations <- function(data) {
  sig_genes <- get_all_genes(data)
  filtered_data <- subset(data, Gene_Name %in% sig_genes)
  cols <- (colnames(filtered_data) == "Gene_Name") | grepl("^Z_", colnames(filtered_data))
  data_wide <- filtered_data[, cols]
  colnames(data_wide) <- gsub("^Z_", "", colnames(data_wide))
  return(data_wide)
}

# ----------------------------
# Heatmap Plotting Function (Single Panel)
# ----------------------------
genes_heatmap_fn <- function(data, plotname) {
  dat <- get_gene_associations(data)
  
  heatmap_data <- dat %>%
    pivot_longer(cols = -Gene_Name, names_to = "Trait", values_to = "Association")
  
  heatmap_data$Gene_Name <- factor(heatmap_data$Gene_Name, levels = rev(sort(unique(heatmap_data$Gene_Name))))
  heatmap_data$Trait <- factor(heatmap_data$Trait, levels = c(
    "HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi", "LVMi", "meanWT",
    "strain_circ", "strain_long", "strain_rad", "NTproBNP", "Troponin"
  ))
  
  p <- ggplot(heatmap_data, aes(x = Trait, y = Gene_Name, fill = Association)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    labs(x = "", y = "", fill = "Association Strength (Z)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"))
  
  ggsave(plotname, p, width = 6.56, height = 12)
  print(p)
}

# ----------------------------
# Plot and Save Single-Panel Heatmap
# ----------------------------
genes_heatmap_fn(mash_post_annotated,
                 "/Users/kbidding/Documents/mashHF/Figures/mash_post_heatmap.png")