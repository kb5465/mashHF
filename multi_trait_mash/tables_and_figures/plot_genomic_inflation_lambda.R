#!/usr/bin/env Rscript

###############################################################
# mashHF: Genomic Inflation Visualization for mashr Results
# Author: Kiran Biddinger
# Description:
# - Computes λ₅₀ (median), λ₉₅, and λ₀₅ genomic inflation stats
# - Plots λ₅₀ as barplots with optional QQ plots
# - Saves figures for mash_st and mash_post datasets
###############################################################

# ----------------------------
# Load Libraries
# ----------------------------
library(ggplot2)
library(dplyr)
library(data.table)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(gaston)

# ----------------------------
# Define Trait Order (for consistent plotting)
# ----------------------------
trait_levels <- c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                  "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                  "NTproBNP", "Troponin")

# ----------------------------
# Function to Calculate λ₅₀ (Median)
# ----------------------------
get_lambda <- function(data) {
  results <- data.frame()
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    z <- data[[paste0("Z_", pheno)]]
    pval <- 2 * pnorm(-abs(z))
    chisq <- z^2
    print(qqplot.pvalues(pval, main = pheno, ylim = c(0, 10)))  # Optional
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    results <- rbind(results, data.frame(Trait = pheno, Inflation = lambda_50))
  }
  results$Trait <- factor(results$Trait, levels = trait_levels)
  return(results)
}

# ----------------------------
# λ₅₀ Barplot Function
# ----------------------------
barplot_fn <- function(data, plotname, title) {
  results <- get_lambda(data)
  mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(length(trait_levels))
  
  p <- ggplot(results, aes(x = Trait, y = Inflation)) +
    geom_hline(yintercept = 1, color = "grey") + 
    geom_bar(aes(fill = Trait), stat = "identity", color = "black", width = 0.9) +
    labs(x = "", y = "Genomic Inflation", title = title) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none") +
    scale_fill_manual(values = mycolors) +
    ylim(0, 2.3)
  
  ggsave(plotname, p, width = 8, height = 5)
  print(p)
}

# ----------------------------
# λ₉₅ Calculator
# ----------------------------
get_lambda_95 <- function(data) {
  results <- data.frame()
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    z <- data[[paste0("Z_", pheno)]]
    pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
    chisq <- qchisq(pval, df = 1, lower.tail = FALSE)
    lambda_95 <- quantile(chisq, 0.95, na.rm = TRUE) / qchisq(0.05, df = 1, lower.tail = FALSE)
    results <- rbind(results, data.frame(Trait = pheno, Inflation = lambda_95))
  }
  results$Trait <- factor(results$Trait, levels = trait_levels)
  return(results)
}

# ----------------------------
# λ₀₅ Calculator
# ----------------------------
get_lambda_5 <- function(data) {
  results <- data.frame()
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    z <- data[[paste0("Z_", pheno)]]
    pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
    chisq <- qchisq(pval, df = 1, lower.tail = FALSE)
    lambda_5 <- quantile(chisq, 0.05, na.rm = TRUE) / qchisq(0.95, df = 1, lower.tail = FALSE)
    results <- rbind(results, data.frame(Trait = pheno, Inflation = lambda_5))
  }
  results$Trait <- factor(results$Trait, levels = trait_levels)
  return(results)
}

# ----------------------------
# Run and Save Plots
# ----------------------------
barplot_fn(mash_st_annotated, 
           "/Users/kbidding/Documents/mashHF/Figures/mash_st_lambda.png", 
           "Single-Trait (λ₅₀)")

barplot_fn(mash_post_annotated, 
           "/Users/kbidding/Documents/mashHF/Figures/mash_post_lambda.png", 
           "Multi-Trait (λ₅₀)")
