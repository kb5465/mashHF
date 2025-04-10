library(gtexr)
library(tidyverse)
library(stats)

# Define the gene names and GENCODE IDs
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

gencode_ids <- c(
  "ENSG00000155657.29", "ENSG00000134571.12", "ENSG00000006747.15", "ENSG00000007306.15", "ENSG00000041515.16",
  "ENSG00000044459.15", "ENSG00000046889.19", "ENSG00000048052.24", "ENSG00000068781.21", "ENSG00000069248.12",
  "ENSG00000072210.19", "ENSG00000080007.8", "ENSG00000083635.8", "ENSG00000086506.3", "ENSG00000088827.13",
  "ENSG00000099899.15", "ENSG00000100344.11", "ENSG00000100360.15", "ENSG00000100441.10", "ENSG00000101331.17",
  "ENSG00000103150.7", "ENSG00000105122.13", "ENSG00000105136.22", "ENSG00000105671.12", "ENSG00000106328.10",
  "ENSG00000108515.18", "ENSG00000109736.15", "ENSG00000109944.11", "ENSG00000112280.18", "ENSG00000115239.24",
  "ENSG00000115525.18", "ENSG00000115705.22", "ENSG00000120279.7", "ENSG00000123064.13", "ENSG00000129696.13",
  "ENSG00000132463.15", "ENSG00000132635.17", "ENSG00000132781.20", "ENSG00000133048.13", "ENSG00000133103.17",
  "ENSG00000133454.16", "ENSG00000133962.8", "ENSG00000134014.17", "ENSG00000134070.5", "ENSG00000136273.13",
  "ENSG00000136783.10", "ENSG00000137198.10", "ENSG00000138346.15", "ENSG00000139988.10", "ENSG00000140400.18",
  "ENSG00000142082.15", "ENSG00000142192.21", "ENSG00000145949.11", "ENSG00000146282.19", "ENSG00000146350.14",
  "ENSG00000148399.13", "ENSG00000149054.16", "ENSG00000151779.13", "ENSG00000154102.11", "ENSG00000157870.17",
  "ENSG00000160191.18", "ENSG00000161958.11", "ENSG00000164404.9", "ENSG00000164483.17", "ENSG00000165071.15",
  "ENSG00000165115.15", "ENSG00000166529.16", "ENSG00000167515.10", "ENSG00000167733.14", "ENSG00000169903.7",
  "ENSG00000170379.21", "ENSG00000170445.16", "ENSG00000171564.12", "ENSG00000176485.13", "ENSG00000178665.16",
  "ENSG00000179044.16", "ENSG00000179142.2", "ENSG00000179270.7", "ENSG00000179709.8", "ENSG00000181031.16",
  "ENSG00000183844.17", "ENSG00000188959.10", "ENSG00000196188.12", "ENSG00000198483.13", "ENSG00000203780.12",
  "ENSG00000205129.9", "ENSG00000213563.7", "ENSG00000214814.7", "ENSG00000235718.9", "ENSG00000242441.8",
  "ENSG00000248385.8", "ENSG00000248919", "ENSG00000249222.1", #"ENSG00000249319", 
  "ENSG00000253910.3",
  "ENSG00000274391.5", "ENSG00000278558.5", "ENSG00000129151.9", "ENSG00000215568.9", "ENSG00000099957.17"
)

# Load required libraries
library(dplyr)
library(tidyr)
library(stats)

# Create the data frame
df <- data.frame(gene_name = genes, gencodeId = gencode_ids, stringsAsFactors = FALSE)

# eQTLs in all outcome GWAS
eqtls_in_all_gwas <- fread("/Users/kbidding/Documents/mashHF/data/GWAS_subset/rsids_in_all_files.csv")

# Tissue names
tissues <- c("Whole_Blood", "Muscle_Skeletal", "Heart_Left_Ventricle", "Heart_Atrial_Appendage")

# Function to retrieve and process eQTLs for a given tissue
process_tissue_eqtls <- function(tissue) {
  # Extract significant eQTLs for the tissue
  eqtls <- get_significant_single_tissue_eqtls(
    gencodeIds = df$gencodeId, 
    tissueSiteDetailIds = tissue, 
    datasetId = "gtex_v10",
    itemsPerPage = 100000
  )
  
  # Remove unnamed SNPs
  eqtls <- subset(eqtls, !is.na(snpId))
  
  # Remove SNPs not in all files
  eqtls <- subset(eqtls, eqtls$snpId %in% eqtls_in_all_gwas$rsid)
  
  # Extract chromosome and position
  eqtls <- separate(eqtls, variantId, into = c("chr", "pos", "a1", "a2", "build"), sep = "_")
  eqtls$chr <- as.numeric(gsub("chr", "", eqtls$chr))
  eqtls$pos <- as.numeric(eqtls$pos)
  
  # Function to remove duplicates within 500000 bp
  remove_duplicates_within_500000_gene <- function(df, gene) {
    df_gene <- subset(df, geneSymbol == gene)
    if (nrow(df_gene) <= 1) return(df_gene)
    
    df_gene <- df_gene[order(df_gene$chr, df_gene$pos), ]
    
    i <- 1
    while (i < nrow(df_gene)) {
      if (df_gene$chr[i] == df_gene$chr[i + 1]) {
        genpos_diff <- abs(df_gene$pos[i + 1] - df_gene$pos[i])
        
        if (genpos_diff < 500000) {
          if (df_gene$pValue[i + 1] > df_gene$pValue[i]) {
            df_gene <- df_gene[-(i + 1), ]  # Remove the second row
          } else {
            df_gene <- df_gene[-i, ]  # Remove the first row
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
  
  # Function to apply filtering to all genes
  remove_duplicates_within_500000 <- function(df) {
    gene_eqtls <- lapply(unique(df$geneSymbol), function(gene) {
      remove_duplicates_within_500000_gene(df, gene)
    })
    return(do.call(rbind, gene_eqtls))
  }
  
  # Apply filtering function
  eqtls_filtered <- remove_duplicates_within_500000(eqtls)
  
  # Rename NES column to beta
  colnames(eqtls_filtered)[colnames(eqtls_filtered) == "nes"] <- "beta"
  
  # Function to compute standard error using lower.tail = TRUE
  compute_standard_error <- function(beta, pValue) {
    if (is.na(beta) | is.na(pValue) | pValue <= 0 | pValue > 1) {
      return(NA)  # Handle invalid cases
    }
    
    adjusted_pValue <- pmax(pValue, 1e-300)
    Z <- -qnorm(adjusted_pValue / 2, lower.tail = TRUE)
    SE <- abs(beta / Z)
    
    return(SE)
  }
  
  # Compute standard errors
  eqtls_filtered$se <- mapply(compute_standard_error, eqtls_filtered$beta, eqtls_filtered$pValue)
  
  # Keep only relevant columns
  eqtls_filtered <- eqtls_filtered[, c("snpId", "geneSymbol", "chr", "pos", "a1", "a2", "beta", "se")]
  
  # Add tissue information
  eqtls_filtered$tissueSiteDetailId <- tissue
  
  return(eqtls_filtered)
}

# Process eQTLs for all tissues and concatenate results
eqtls_all_tissues <- do.call(rbind, lapply(tissues, process_tissue_eqtls))

# Final concatenated dataframe with all tissue results
eqtls_all_tissues

# Write output
write.csv(eqtls_all_tissues, "/Users/kbidding/Documents/mashHF/data/mashHF_eqtl_031825.csv")

