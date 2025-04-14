#!/usr/bin/env Rscript

###############################################################
# mashHF: Extract N_CARRIERS by Subset for Canonical LOF Genes
# Author: Kiran Biddinger
# Description:
# - Reads canonical LOF result files for multiple phenotypes
# - Extracts N_CARRIERS per gene per subset
# - Annotates genes with HGNC names
# - Filters to a predefined gene set of interest
###############################################################

# ----------------------------
# 1. Load Required Libraries
# ----------------------------
library(tidyverse)
library(data.table)
library(biomaRt)

# ----------------------------
# 2. Define Input Paths and File List
# ----------------------------
data_dir <- "/Users/kbidding/Documents/mashHF/data/st"
files <- list.files(data_dir, pattern = "_canonical_LOF_0.01_merged.tsv$", full.names = TRUE)

# Extract phenotype names (e.g., "HF", "LVEF", "NTproBNP")
phenos <- gsub("_canonical_LOF_0.01_merged.tsv", "", basename(files))
names(files) <- phenos

# ----------------------------
# 3. Load Files and Extract N_CARRIERS
# ----------------------------
read_and_compute_ncarriers <- function(file, pheno) {
  df <- as.data.frame(fread(file))
  df[[pheno]] <- df$N_CARRIERS  # Assign to new column with phenotype name
  df_subset <- df[, c("ID", pheno)]
  return(df_subset)
}

ncar_list <- map2(files, phenos, read_and_compute_ncarriers)
ncar_data <- reduce(ncar_list, full_join, by = "ID")  # Merge on gene ID

# ----------------------------
# 4. Annotate Ensembl Gene IDs with HGNC Symbols
# ----------------------------
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ncar_data$ID,
  mart = ensembl
)

colnames(gene_annotations) <- c("ID", "Gene_Name")
ncar_data <- left_join(ncar_data, gene_annotations, by = "ID")

# ----------------------------
# 5. Filter to Genes of Interest
# ----------------------------
genes <- c(
  "TTN", "MYBPC3", "SCIN", "CEACAM7", "MYO16", "CNTLN", "PREX2", "HDAC9",
  "STON1-GTF2A1L", "NUP133", "ALDH3A2", "DDX43", "NUFIP1", "HBQ1", "SIGLEC1", "TRMT2A",
  "PNPLA3", "IFT27", "KHNYN", "CCM2L", "MLYCD", "RASAL3", "ZNF419", "DDX49",
  "FSCN3", "ENO3", "MFSD10", "JHY", "COL9A1", "ASB3", "ST3GAL5", "TPO",
  "MYCT1", "DDX54", "TTI2", "GRSF1", "PCED1A", "MUTYH", "CHI3L1", "COG6",
  "MYO18B", "CATSPERB", "ELP3", "IRAK2", "HUS1", "NIPSNAP3A", "GMPR", "DNA2",
  "RDH12", "MAN2C1", "SIRT3", "APP", "MYLK4", "RARS2", "TBC1D32", "DPH7",
  "ZNF215", "NBAS", "C16orf74", "PRXL2B", "PDE9A", "FGF11", "GDF9", "SAMD3",
  "TMEM71", "KIF27", "ZSCAN21", "TRAPPC2L", "HSD11B1L", "TM4SF4", "TCAF2", "HARS1",
  "FGB", "PLAAT3", "ZNF713", "EXOC3L1", "CYP11B2", "PCARE", "NLRP8", "RPH3AL",
  "FAM3B", "C9orf152", "CTSE", "ANKRD35", "FANK1", "CFAP96", "C8orf82", "FER1L6",
  "MFRP", "GTF2A1L", "TARM1", "ATP5MF-PTCD1", "ATP5MGL", "PCDHGB2", "TPTE", "TMEM191B",
  "BBOX1", "GAB4", "P2RX6"
)

ncar_data <- filter(ncar_data, Gene_Name %in% genes)

# ----------------------------
# 6. Format Final Table
# ----------------------------
# Rename columns for clarity
ncar_final <- data.frame(
  Gene = ncar_data$Gene_Name,
  All = round(ncar_data$HF),
  LV_Subset = round(ncar_data$LVEF),
  LV_ML_Subset = round(ncar_data$LVMi),
  Proteomics_Subset = round(ncar_data$NTproBNP)
)

# Optional: write to file
# write.csv(ncar_final, "N_carriers_by_subset.csv", row.names = FALSE)

# View output
print(head(ncar_final))