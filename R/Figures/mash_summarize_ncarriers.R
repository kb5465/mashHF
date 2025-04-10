# ---- Load Libraries ----
library(tidyverse)
library(data.table)
library(biomaRt)

# ---- File Settings ----
# Path to directory with TSVs
data_dir <- "/Users/kbidding/Documents/mashHF/data/st"

# Get list of all canonical_LOF result files
files <- list.files(data_dir, pattern = "_canonical_LOF_0.01_merged.tsv$", full.names = TRUE)

# Extract phenotype names (PHENO)
phenos <- gsub("_canonical_LOF_0.01_merged.tsv", "", basename(files))
names(files) <- phenos

# ---- Load All Files and Compute N_carriers ----
# Function to read file and compute z-scores
read_and_compute_ncarriers <- function(file, pheno) {
  df <- as.data.frame(fread(file))
  
  # Compute Z-scores as beta / se
  df[[pheno]] <- df$N_CARRIERS
  
  # Keep only gene ID and Z-score
  df_subset <- df[, c("ID", pheno)]
  return(df_subset)
}

# Apply the function to all files
ncar_list <- map2(files, phenos, read_and_compute_ncarriers)

# Merge all data frames by gene ID
ncar_data <- reduce(ncar_list, full_join, by = "ID")

# ---- Annotate Genes with Gene Names ----
# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch gene annotations
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ncar_data$ID,
  mart = ensembl
)

# Rename for merging
colnames(gene_annotations) <- c("ID", "Gene_Name")

# Merge to get gene names
ncar_data <- left_join(ncar_data, gene_annotations, by = "ID")


# Define relevant genes
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

# Select relevant genes
ncar_data <- subset(ncar_data, Gene_Name %in% genes)

# Trim to unique population-specific values
ncar_final <- data.frame(Gene = ncar_data$Gene_Name,
                         All = round(ncar_data$HF),
                         LV_Subset = round(ncar_data$LVEF),
                         LV_ML_Subset = round(ncar_data$LVMi),
                         Proteomics_Subset = round(ncar_data$NTproBNP))






