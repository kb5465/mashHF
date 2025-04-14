library(tidyverse)
library(data.table)
library(stats)

# Load eqtl data
eqtls <- fread("/Users/kbidding/Documents/mashHF/data/mashHF_eqtl_031825.csv")
genes <- unique(eqtls$geneSymbol)
eqtl_snps <- unique(eqtls$snpId)

# Load GWAS files
gwas_file_names <- list.files("/Users/kbidding/Documents/mashHF/data/GWAS_subset")
phenos <- gsub(".csv", "", gwas_file_names)

# Tissue types
tissues <- c("Whole_Blood", "Muscle_Skeletal", "Heart_Left_Ventricle", "Heart_Atrial_Appendage")

# Perform Mendelian randomization for each gene and each outcome, within each tissue
outcome_associated_eqtls <- function(gene, outcome, tissue) {
  # Filter to gene, tissue
  exposure_data <- subset(eqtls, geneSymbol == gene & tissueSiteDetailId == tissue)
  exposure_data <- exposure_data[, c("snpId", "a1", "a2", "beta", "se", "geneSymbol")]
  colnames(exposure_data)[2:5] <- paste0(colnames(exposure_data)[2:5], "_exposure")

  # Get relevant outcome data
  outcome_data <- fread(paste0("/Users/kbidding/Documents/mashHF/data/GWAS_subset/", outcome, ".csv"))
  outcome_data$V1 <- NULL
  colnames(outcome_data)[2:5] <- paste0(colnames(outcome_data)[2:5], "_outcome")
  z <- outcome_data$beta_outcome / outcome_data$se_outcome
  outcome_data$pval_outcome <- 2 * pnorm(-abs(z))
  
  # Merge files on SNP ID
  merged_data <- merge(exposure_data, outcome_data, by.x = "snpId", by.y = "rsid")
  if(nrow(merged_data) == 0) {return()}
  
  # Clean merged data
  merged_data <- subset(merged_data, (a1_exposure == ea_outcome & a2_exposure == nea_outcome) | 
                      (a2_exposure == ea_outcome & a1_exposure == nea_outcome))
  
  # Align flipped alleles
  merged_data$beta_outcome[merged_data$a1_exposure == merged_data$nea_outcome] <- -1 * merged_data$beta_outcome[merged_data$a1_exposure == merged_data$nea_outcome]
  merged_data$ea_outcome[merged_data$a1_exposure == merged_data$nea_outcome] <- merged_data$a2_exposure[merged_data$a1_exposure == merged_data$nea_outcome]
  merged_data$ea_outcome[merged_data$a1_exposure == merged_data$nea_outcome] <- merged_data$a1_exposure[merged_data$a1_exposure == merged_data$nea_outcome]
  
  # Print significant associations
  if(any(merged_data$pval_outcome < 0.05 / 85)) {
    #print(sprintf("Gene: %s, Outcome: %s, Tissue: %s", gene, outcome, tissue))
    #print(merged_data)
  }
  
  # Calculate tissue-level cauchy combination p-value
  cauchy_pval <- CCT(merged_data$pval_outcome)
  
  # Return data
  return(cauchy_pval)
}

# Find outcome associated eqtls
for (gene in genes) {
  for (pheno in phenos) {
    cauchy_tissue_pvals <- c()
    for (tissue in tissues) {
      cauchy_tissue_pvals <- c(cauchy_tissue_pvals, outcome_associated_eqtls(gene, pheno, tissue))
    }
    if (length(cauchy_tissue_pvals) > 0) {
      cauchy_trait_pval <- CCT(cauchy_tissue_pvals)
      if (cauchy_trait_pval < 0.05) {
        print(sprintf("Gene: %s, Outcome: %s: Cauchy p-value = %f", gene, pheno, cauchy_trait_pval))
      }
    }
  }
}


# From: https://rdrr.io/github/xihaoli/STAAR/src/R/CCT.R
CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }
  
  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }
  
  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }
  
  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }
  
  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}


