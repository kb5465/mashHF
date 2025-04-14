#!/usr/bin/env Rscript

###############################################################
# mashHF: Count Bonferroni-significant Gene Hits per Trait
# Author: Kiran Biddinger
# Description:
# - Computes number of genes with p < 0.05 / m for each trait
# - Plots barplots for mash_prior and mash_post data
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
# Define Trait Order (for consistency)
# ----------------------------
trait_order <- c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                 "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                 "NTproBNP", "Troponin")

# ----------------------------
# Function to Count Significant Genes per Trait
# ----------------------------
get_nhits <- function(data) {
  results <- data.frame()
  
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    z <- data[[paste0("Z_", pheno)]]
    pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
    bonf <- 0.05 / sum(!is.na(pval))
    nhits <- sum(pval < bonf)
    
    results <- rbind(results, data.frame(Trait = pheno, Inflation = nhits))
  }
  
  results$Trait <- factor(results$Trait, levels = trait_order)
  return(results)
}

# ----------------------------
# Function to Plot Number of Hits per Trait
# ----------------------------
barplot_nhits_fn <- function(data, plotname, title) {
  results <- get_nhits(data)
  mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(length(trait_order))
  
  p <- ggplot(results, aes(x = Trait, y = Inflation)) +
    geom_hline(yintercept = 0, color = "grey") + 
    geom_bar(aes(fill = Trait), stat = "identity", color = "black", width = 0.9) +
    labs(x = "", y = "Number of Genes Associated", title = title) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1), 
          axis.title.x = element_text(size = 11, color = "black"),  
          axis.text.y = element_text(size = 11, color = "black"), 
          axis.title.y = element_text(size = 11, color = "black"),
          axis.ticks.y = element_line(color = "black"),
          axis.ticks.x = element_line(color = "black"),
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_fill_manual(values = mycolors)
  
  ggsave(plotname, p, width = 8, height = 5)
  print(p)
}

# ----------------------------
# Generate and Save Plots
# ----------------------------
barplot_nhits_fn(mash_prior_annotated, 
                 "/Users/kbidding/Documents/mashHF/Figures/mash_prior_nhits.png", 
                 "Single-Trait")

barplot_nhits_fn(mash_post_annotated, 
                 "/Users/kbidding/Documents/mashHF/Figures/mash_post_nhits.png", 
                 "Multi-Trait")