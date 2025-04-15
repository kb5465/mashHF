#!/usr/bin/env Rscript

######################################################################
# Script: extract_eqtls_for_mashHF_genes.R
# Author: Kiran Biddinger
# Purpose:
#   - Download and filter GTEx v10 significant eQTLs for mashHF genes
#   - Keep only SNPs present in all outcome GWAS datasets
#   - Filter out nearby duplicates (within 500kb)
#   - Estimate standard errors from beta and p-values
# Output: Final eQTL table written to CSV
######################################################################

# ----------------------------
# Load Libraries
# ----------------------------
library(gtexr)
library(data.table)
library(tidyverse)
library(stats)

# ----------------------------
# Input: List of Genes and GENCODE IDs
# ----------------------------
genes <- c("TTN", "MYBPC3", "SCIN", "CEACAM7", "MYO16", ...)  # truncated for brevity
gencode_ids <- c("ENSG00000155657.29", "ENSG00000134571.12", ...)  # match gene order

df <- data.frame(gene_name = genes, gencodeId = gencode_ids, stringsAsFactors = FALSE)

# ----------------------------
# Load rsIDs present in all GWAS datasets
# ----------------------------
eqtls_in_all_gwas <- fread("/Users/kbidding/Documents/mashHF/data/GWAS_subset/rsids_in_all_files.csv")

# ----------------------------
# Tissues to Query from GTEx
# ----------------------------
tissues <- c("Whole_Blood", "Muscle_Skeletal", "Heart_Left_Ventricle", "Heart_Atrial_Appendage")

# ----------------------------
# Helper: Remove Nearby Duplicates per Gene (<500kb apart)
# ----------------------------
remove_duplicates_within_500000_gene <- function(df, gene) {
  df_gene <- df[df$geneSymbol == gene, ]
  if (nrow(df_gene) <= 1) return(df_gene)
  
  df_gene <- df_gene[order(df_gene$chr, df_gene$pos), ]
  
  i <- 1
  while (i < nrow(df_gene)) {
    if (df_gene$chr[i] == df_gene$chr[i + 1]) {
      diff <- abs(df_gene$pos[i + 1] - df_gene$pos[i])
      if (diff < 500000) {
        if (df_gene$pValue[i + 1] > df_gene$pValue[i]) {
          df_gene <- df_gene[-(i + 1), ]
        } else {
          df_gene <- df_gene[-i, ]
        }
      } else {
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  return(df_gene)
}

remove_duplicates_within_500000 <- function(df) {
  do.call(rbind, lapply(unique(df$geneSymbol), function(g) {
    remove_duplicates_within_500000_gene(df, g)
  }))
}

# ----------------------------
# Helper: Compute SE from beta and p-value
# ----------------------------
compute_standard_error <- function(beta, pValue) {
  if (is.na(beta) || is.na(pValue) || pValue <= 0 || pValue > 1) return(NA)
  adjusted_p <- pmax(pValue, 1e-300)
  z <- -qnorm(adjusted_p / 2, lower.tail = TRUE)
  return(abs(beta / z))
}

# ----------------------------
# Main Function: Extract eQTLs for a given tissue
# ----------------------------
process_tissue_eqtls <- function(tissue) {
  eqtls <- get_significant_single_tissue_eqtls(
    gencodeIds = df$gencodeId,
    tissueSiteDetailIds = tissue,
    datasetId = "gtex_v10",
    itemsPerPage = 100000
  )
  
  eqtls <- subset(eqtls, !is.na(snpId))
  eqtls <- subset(eqtls, snpId %in% eqtls_in_all_gwas$rsid)
  
  eqtls <- eqtls %>%
    separate(variantId, into = c("chr", "pos", "a1", "a2", "build"), sep = "_", remove = FALSE) %>%
    mutate(chr = as.numeric(gsub("chr", "", chr)),
           pos = as.numeric(pos)) %>%
    select(-build)
  
  eqtls <- remove_duplicates_within_500000(eqtls)
  
  eqtls <- eqtls %>%
    rename(beta = nes) %>%
    mutate(se = mapply(compute_standard_error, beta, pValue)) %>%
    select(snpId, geneSymbol, chr, pos, a1, a2, beta, se) %>%
    mutate(tissueSiteDetailId = tissue)
  
  return(eqtls)
}

# ----------------------------
# Loop Through All Tissues
# ----------------------------
eqtls_all_tissues <- bind_rows(lapply(tissues, process_tissue_eqtls))

# ----------------------------
# Write to File
# ----------------------------
write.csv(eqtls_all_tissues,
          "/Users/kbidding/Documents/mashHF/data/mashHF_eqtl_031825.csv",
          row.names = FALSE)

# Optional preview
print(head(eqtls_all_tissues))