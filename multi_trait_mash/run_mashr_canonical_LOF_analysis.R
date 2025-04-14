#!/usr/bin/env Rscript

###############################################
# mashHF: mashr Model on Canonical LOF 0.01 Variants
# Author: Kiran Biddinger
# Description:
# - Load gene-level RVAT stats (BETA, SE)
# - Fit mashr model using canonical + data-driven covariances
# - Identify significant hits pre- and post-mash
# - Annotate and export results
###############################################

# Load required packages
library(data.table)
library(mashr)
library(tidyverse)
library(biomaRt)
library(openxlsx)

# ----------------------------
# 1. Load and Filter Input Data
# ----------------------------
rvat <- fread("/Users/kbidding/Documents/mashHF/data/mashHF_canonical_LOF_0.01_merged_firth.tsv")
mash_st <- drop_na(rvat) |> as.data.frame()

# ----------------------------
# 2. Prepare mashr Input
# ----------------------------
Bhat <- as.matrix(mash_st[, grep("^BETA_", colnames(mash_st))])
Shat <- as.matrix(mash_st[, grep("^SE_", colnames(mash_st))])
mash_data <- mash_set_data(Bhat, Shat)

# Compute Z-scores
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  mash_st[[paste0("Z_", pheno)]] <- mash_st[[paste0("BETA_", pheno)]] / mash_st[[paste0("SE_", pheno)]]
}

# ----------------------------
# 3. Identify Strong Tests
# ----------------------------
m_1by1 <- mash_1by1(mash_data)
Vhat <- estimate_null_correlation_simple(mash_data)
strong_subset <- get_significant_results(m_1by1, 0.1)
mash_strong <- mash_set_data(Bhat[strong_subset,], Shat[strong_subset,], V = Vhat)

# ----------------------------
# 4. Define Covariance Matrices
# ----------------------------
U_pca <- cov_pca(mash_strong, 5)
U_ed  <- cov_ed(mash_strong, U_pca)
U_cov <- cov_canonical(mash_data)

# ----------------------------
# 5. Run mashr Model
# ----------------------------
mash_result <- mash(mash_data, Ulist = c(U_ed, U_cov), pi_thresh = 1e-30)

# ----------------------------
# 6. Extract Posterior Estimates
# ----------------------------
mash_post <- as.data.frame(mash_result$result$PosteriorMean)
mash_post_sd <- as.data.frame(mash_result$result$PosteriorSD)
mash_post_lfsr <- as.data.frame(mash_result$result$lfsr)

# Add IDs
mash_post$ID <- mash_st$ID
mash_post_sd$ID <- mash_st$ID
mash_post_lfsr$ID <- mash_st$ID

# Rename columns
colnames(mash_post_sd)   <- gsub("BETA", "SE", colnames(mash_post_sd))
colnames(mash_post_lfsr) <- gsub("BETA", "lfsr", colnames(mash_post_lfsr))

# Merge all posteriors
mash_post <- mash_post |>
  inner_join(mash_post_sd, by = "ID") |>
  inner_join(mash_post_lfsr, by = "ID")

# ----------------------------
# 7. Compute Z-scores (Post-Mash)
# ----------------------------
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  mash_post[[paste0("Z_", pheno)]] <- mash_post[[paste0("BETA_", pheno)]] / mash_post[[paste0("SE_", pheno)]]
}

# ----------------------------
# 8. Gene Annotation
# ----------------------------
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = mash_post$ID,
  mart = ensembl
) |> rename(ID = ensembl_gene_id, Gene_Name = external_gene_name)

# Merge annotations
mash_post_annotated  <- left_join(mash_post, gene_annotations, by = "ID")
mash_st_annotated <- left_join(mash_st, gene_annotations, by = "ID")

# ----------------------------
# 9. Identify Significant Hits (Bonferroni)
# ----------------------------
bonf_sig_thr <- 0.05 / nrow(mash_st)
z_thr <- qnorm(1 - bonf_sig_thr / 2)

# Export single-trait-based significant hits
wb_st <- createWorkbook()
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  sig_hits <- mash_st_annotated |> filter(abs(.data[[paste0("Z_", pheno)]]) > z_thr)
  sig_hits$PVAL <- 2 * (1 - pnorm(abs(sig_hits[[paste0("Z_", pheno)]])))
  sig_hits <- sig_hits[, c("Gene_Name", paste0("BETA_", pheno), paste0("SE_", pheno), "PVAL")]
  addWorksheet(wb_st, pheno)
  writeData(wb_st, pheno, sig_hits)
}
saveWorkbook(wb_st, "/Users/kbidding/Documents/mashHF/Tables/mash_st_significant_hits.xlsx", overwrite = TRUE)

# Export mash-based significant hits
wb_post <- createWorkbook()
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  sig_hits <- mash_post_annotated |> filter(abs(.data[[paste0("Z_", pheno)]]) > z_thr)
  sig_hits$PVAL <- 2 * (1 - pnorm(abs(sig_hits[[paste0("Z_", pheno)]])))
  sig_hits <- sig_hits[, c("Gene_Name", paste0("BETA_", pheno), paste0("SE_", pheno), "PVAL", paste0("lfsr_", pheno))]
  addWorksheet(wb_post, pheno)
  writeData(wb_post, pheno, sig_hits)
}
saveWorkbook(wb_post, "/Users/kbidding/Documents/mashHF/Tables/mash_post_significant_hits.xlsx", overwrite = TRUE)

# ----------------------------
# 10. Gene-by-Gene Comparison for Key Cardiomyopathy Genes
# ----------------------------
genes_highlight <- c("TTN", "BAG3", "DES", "FLNC", "LMNA", "MYH7", "PLN", "RBM20", "SCN5A", "TNNC1", "TNNT2", "DSP", "MYBPC3")

for (gene in genes_highlight) {
  st_row <- mash_st_annotated |> filter(Gene_Name == gene)
  post_row  <- mash_post_annotated  |> filter(Gene_Name == gene)
  
  if (nrow(st_row) == 0 | nrow(post_row) == 0) {
    message("Skipping gene: ", gene, " (not found)")
    next
  }
  
  common_cols <- intersect(
    colnames(st_row)[grepl("BETA|SE|Z", colnames(st_row))],
    colnames(post_row)[grepl("BETA|SE|Z", colnames(post_row))]
  )
  
  combined <- rbind(st_row[1, common_cols], post_row[1, common_cols]) |> t() |> as.data.frame()
  colnames(combined) <- c("Naive (single-trait)", "Posterior (multi-trait)")
  combined$Metric <- rownames(combined)
  
  cat("\n--- Comparison for gene:", gene, "---\n")
  print(combined)
}
