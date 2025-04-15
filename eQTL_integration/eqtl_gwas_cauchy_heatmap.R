#!/usr/bin/env Rscript

################################################################################
# Title: Pan-Tissue eQTL × GWAS Association Mapping with Cauchy Combination
# Description:
#   Combines GTEx cis-eQTL data with GWAS summary statistics for cardiac traits.
#   For each gene × trait × tissue, computes a Cauchy combined p-value. Then
#   aggregates p-values across tissues. Significant associations are visualized.
#
# Author: Kiran Biddinger
# Date: April 2025
################################################################################

# ---- Load Libraries ----
library(tidyverse)
library(data.table)
library(stats)

# ---- Load and Preprocess eQTL Data ----
eqtls <- fread("~/Documents/mashHF/data/mashHF_eqtl_031825.csv")
eqtls <- subset(eqtls, geneSymbol != "ANKRD35")  # Drop gene with NA info

genes <- unique(eqtls$geneSymbol)
eqtl_snps <- unique(eqtls$snpId)

# ---- GWAS File List ----
gwas_dir <- "~/Documents/mashHF/data/GWAS_subset"
gwas_file_names <- list.files(gwas_dir, pattern = "\\.csv$", full.names = TRUE)
gwas_file_names <- gwas_file_names[!grepl("rsids_in_all_files", gwas_file_names)]
phenos <- gsub("\\.csv$", "", basename(gwas_file_names))

# ---- Tissues ----
tissues <- c("Whole_Blood", "Muscle_Skeletal", "Heart_Left_Ventricle", "Heart_Atrial_Appendage")

# ---- Define Cauchy Combination Test ----
CCT <- function(pvals, weights = NULL) {
  pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0) return(NA)
  if (any(pvals == 0)) return(0)
  if (any(pvals == 1)) return(1)
  if (is.null(weights)) weights <- rep(1 / length(pvals), length(pvals))
  weights <- weights / sum(weights)
  is.small <- pvals < 1e-16
  stat <- if (any(is.small)) {
    sum((weights[is.small] / pvals[is.small]) / pi) + 
      sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  } else {
    sum(weights * tan((0.5 - pvals) * pi))
  }
  if (stat > 1e15) return((1 / stat) / pi)
  return(1 - pcauchy(stat))
}

# ---- Function: Gene × Trait × Tissue Association ----
outcome_associated_eqtls <- function(gene, outcome, tissue) {
  exposure_data <- eqtls %>%
    filter(geneSymbol == gene, tissueSiteDetailId == tissue) %>%
    select(snpId, a1, a2, beta, se) %>%
    rename_with(~ paste0(., "_exposure"), c("a1", "a2", "beta", "se"))
  
  if (nrow(exposure_data) == 0) return(NA)
  
  outcome_data <- fread(file.path(gwas_dir, paste0(outcome, ".csv"))) %>%
    select(-V1) %>%
    rename_with(~ paste0(., "_outcome"), c("ea", "nea", "beta", "se")) %>%
    mutate(pval_outcome = 2 * pnorm(-abs(beta_outcome / se_outcome)))
  
  merged <- merge(exposure_data, outcome_data, by.x = "snpId", by.y = "rsid")
  if (nrow(merged) == 0) return(NA)
  
  merged <- merged %>% 
    filter((a1_exposure == ea_outcome & a2_exposure == nea_outcome) |
             (a2_exposure == ea_outcome & a1_exposure == nea_outcome))
  
  merged$beta_outcome[merged$a1_exposure == merged$nea_outcome] <-
    -merged$beta_outcome[merged$a1_exposure == merged$nea_outcome]
  
  if (nrow(merged) == 0) return(NA)
  
  return(CCT(merged$pval_outcome))
}

# ---- Per-Tissue p-values ----
tissue_results <- list()
for (t in tissues) {
  message("Processing tissue: ", t)
  tissue_df <- data.frame(gene = genes)
  for (p in phenos) {
    message(sprintf(" - Trait: %s", p))
    tissue_df[[p]] <- sapply(genes, function(g) outcome_associated_eqtls(g, p, t))
  }
  tissue_results[[t]] <- tissue_df
}

# ---- Pan-Tissue Aggregation ----
final_df <- data.frame(gene = genes)

# Add eQTL counts
final_df$num_eqtls_all_tissues <- sapply(genes, function(g) {
  sum(sapply(tissues, function(t) nrow(filter(eqtls, geneSymbol == g, tissueSiteDetailId == t))))
})
final_df$num_eqtls_all_tissues_norep <- sapply(genes, function(g) {
  length(unique(eqtls$snpId[eqtls$geneSymbol == g]))
})

# Aggregate across tissues using CCT
for (p in phenos) {
  message("Aggregating across tissues for trait: ", p)
  final_df[[p]] <- sapply(genes, function(g) {
    pvals <- sapply(tissues, function(t) {
      df <- tissue_results[[t]]
      val <- df[df$gene == g, p]
      if (length(val) == 0) return(NA)
      return(val)
    })
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) == 0) return(NA)
    return(CCT(pvals))
  })
}

# ---- Identify Genes with Significant Trait Associations ----
sig_final_df <- final_df %>%
  filter(apply(select(., all_of(phenos)), 1, function(x) sum(x < 0.05 / 85, na.rm = TRUE) > 0)) %>%
  select(-num_eqtls_all_tissues, -num_eqtls_all_tissues_norep) %>%
  mutate(across(HFpEF:HFrEF, ~ NULL))  # Remove unwanted traits if needed

# ---- Prepare Heatmap Data ----
long_df <- sig_final_df %>%
  pivot_longer(-gene, names_to = "Trait", values_to = "pval") %>%
  filter(!is.na(pval)) %>%
  mutate(
    Association = -log10(pval),
    Gene_Name = factor(gene, levels = sort(unique(gene), decreasing = TRUE)),
    Side = ifelse(Trait %in% c("DCM", "HCM", "HF", "HFpEF", "HFrEF"), "Disease", "Imaging")
  )

# ---- Plot Heatmap ----
ggplot(long_df, aes(x = Trait, y = Gene_Name, fill = Association)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "darksalmon", mid = "thistle", high = "darkblue",
    midpoint = -log10(0.05), name = "-log10(p)"
  ) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  )

# ---- Print Gene–Trait Pairs with FDR-Significant Associations ----
for (g in genes) {
  for (p in phenos) {
    pval <- final_df[final_df$gene == g, p]
    if (!is.na(pval) && pval < (0.05 / 85)) {
      message(sprintf("Gene: %s, Trait: %s, Cauchy p = %g", g, p, pval))
    }
  }
}