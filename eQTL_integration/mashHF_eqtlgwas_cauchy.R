library(tidyverse)
library(data.table)
library(stats)

# Load eqtl data
eqtls <- fread("/Users/kbidding/Documents/mashHF/data/mashHF_eqtl_031825.csv")
genes <- unique(eqtls$geneSymbol)
genes <- subset(genes, genes != "ANKRD35") # NA Values
eqtl_snps <- unique(eqtls$snpId)

# Load GWAS file names and phenotype names
gwas_file_names <- list.files("/Users/kbidding/Documents/mashHF/data/GWAS_subset")
gwas_file_names <- subset(gwas_file_names, gwas_file_names != "rsids_in_all_files.csv")
phenos <- gsub(".csv", "", gwas_file_names)

# Tissue types
tissues <- c("Whole_Blood", "Muscle_Skeletal", "Heart_Left_Ventricle", "Heart_Atrial_Appendage")

# Cauchy combination test
CCT <- function(pvals, weights=NULL){
  pvals <- subset(pvals, !is.na(pvals))
  if(sum(is.na(pvals)) > 0) stop("Cannot have NAs in the p-values!")
  if((sum(pvals < 0) + sum(pvals > 1)) > 0) stop("All p-values must be between 0 and 1!")
  if(any(pvals == 0)) return(0)
  if(any(pvals == 1)) {
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  if(is.null(weights)) {
    weights <- rep(1 / length(pvals), length(pvals))
  } else {
    if(length(weights) != length(pvals)) stop("Weights length must match p-values!")
    weights <- weights / sum(weights)
  }
  is.small <- (pvals < 1e-16)
  cct.stat <- if (sum(is.small) == 0) {
    sum(weights * tan((0.5 - pvals) * pi))
  } else {
    sum((weights[is.small] / pvals[is.small]) / pi) + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  }
  if(cct.stat > 1e+15) {
    return((1 / cct.stat) / pi)
  } else {
    return(1 - pcauchy(cct.stat))
  }
}

# Function to get Cauchy p-value for gene, phenotype, and tissue
outcome_associated_eqtls <- function(gene, outcome, tissue) {
  exposure_data <- subset(eqtls, geneSymbol == gene & tissueSiteDetailId == tissue)
  if (nrow(exposure_data) == 0) return(NA)
  num_eqtls <- nrow(exposure_data)
  exposure_data <- exposure_data[, c("snpId", "a1", "a2", "beta", "se", "geneSymbol")]
  colnames(exposure_data)[2:5] <- paste0(colnames(exposure_data)[2:5], "_exposure")
  
  outcome_data <- fread(paste0("/Users/kbidding/Documents/mashHF/data/GWAS_subset/", outcome, ".csv"))
  outcome_data$V1 <- NULL
  colnames(outcome_data)[2:5] <- paste0(colnames(outcome_data)[2:5], "_outcome")
  z <- outcome_data$beta_outcome / outcome_data$se_outcome
  outcome_data$pval_outcome <- 2 * pnorm(-abs(z))
  
  merged_data <- merge(exposure_data, outcome_data, by.x = "snpId", by.y = "rsid")
  if (nrow(merged_data) == 0) return(NA)
  
  merged_data <- subset(merged_data, (a1_exposure == ea_outcome & a2_exposure == nea_outcome) | 
                          (a2_exposure == ea_outcome & a1_exposure == nea_outcome))
  
  merged_data$beta_outcome[merged_data$a1_exposure == merged_data$nea_outcome] <- 
    -1 * merged_data$beta_outcome[merged_data$a1_exposure == merged_data$nea_outcome]
  
  if (nrow(merged_data) == 0) return(NA)
  
  return(CCT(merged_data$pval_outcome))
}

# Create per-tissue data frames
tissue_results <- list()

for (tissue in tissues) {
  tissue_df <- data.frame(gene = genes, stringsAsFactors = FALSE)
  tissue_df$num_eqtls <- sapply(genes, function(g) {
    nrow(subset(eqtls, geneSymbol == g & tissueSiteDetailId == tissue))
    })
  
  for (pheno in phenos) {
    cat(sprintf("Processing %s - %s\n", tissue, pheno))
    pvals <- sapply(genes, function(g) outcome_associated_eqtls(g, pheno, tissue))
    tissue_df[[pheno]] <- pvals
  }
  
  tissue_results[[tissue]] <- tissue_df
}

# View result for one tissue
# print(tissue_results[["Heart_Left_Ventricle"]])

# Create final summary table
final_df <- data.frame(gene = genes, stringsAsFactors = FALSE)

# Sum of eQTLs across tissues (counting repeats)
final_df$num_eqtls_all_tissues <- sapply(genes, function(g) {
  sum(sapply(tissues, function(t) {
    nrow(subset(eqtls, geneSymbol == g & tissueSiteDetailId == t))
  }))
})

# Sum of eQTLs across tissues (not counting repeats)
final_df$num_eqtls_all_tissues_norep <- sapply(genes, function(g) {
  length(unique(eqtls$snpId[eqtls$geneSymbol == g]))
})

# Combine p-values across tissues for each phenotype using CCT
for (pheno in phenos) {
  cat(sprintf("Computing pan-tissue p-values for phenotype: %s\n", pheno))
  pan_pvals <- sapply(genes, function(g) {
    pvals <- sapply(tissues, function(t) {
      df <- tissue_results[[t]]
      val <- df[df$gene == g, pheno]
      if (length(val) == 0) return(NA)
      return(val)
    })
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) == 0) return(NA)
    return(CCT(pvals))
  })
  final_df[[pheno]] <- pan_pvals
}

# View or write the final table
print(head(final_df))

# Subset to rows with >1 nominal association
sig_final_df <- subset(final_df, sapply(genes, function(g) {
  sum(final_df[final_df$gene == g, 3:ncol(final_df)] < 0.05 / 80, na.rm = T) > 0
  }))
sig_final_df$num_eqtls_all_tissues <- NULL
sig_final_df$num_eqtls_all_tissues_norep <- NULL
sig_final_df$HFpEF <- NULL
sig_final_df$HFrEF <- NULL
# Optional: write to CSV
#write.csv(final_df, "/Users/kbidding/Documents/mashHF/Tables/pan_tissue_summary.csv", row.names = FALSE)

# Reshape to long format
long_df <- sig_final_df %>%
  pivot_longer(
    cols = -gene,
    names_to = "Trait",
    values_to = "pval"
  ) %>%
  filter(!is.na(pval)) %>%
  mutate(
    Association = -log10(pval),
    Gene_Name = gene,
    Side = ifelse(Trait %in% c("DCM", "HCM", "HF", "HFpEF", "HFrEF"), "Disease", "Imaging")
  )

# Factor gene names alphabetically (top to bottom)
long_df <- long_df %>%
  mutate(Gene_Name = factor(Gene_Name, levels = sort(unique(Gene_Name), decreasing = TRUE)))

# Plot heatmap
ggplot(long_df, aes(x = Trait, y = Gene_Name, fill = Association)) +
  geom_tile() +
  scale_fill_gradient2(low = "darksalmon", mid = "thistle", high = "darkblue", midpoint = -log10(0.05),
                       name = "-log10(p)") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18), 
    axis.text.y = element_text(size = 18), 
    axis.title.y = element_blank(),
    strip.text.y.left = element_text(size = 12, angle = 90),
    strip.placement = "outside",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    panel.spacing.x = unit(0.5, "cm"),
    panel.grid = element_blank() 
  )


# Find outcome associated eqtls
for (g in genes) {
  for (p in phenos) {
    cauchy_pval <- final_df[final_df$gene == g, p]
    if (identical(cauchy_pval, numeric(0)) | is.na(cauchy_pval)) {
      next
    }
    if (cauchy_pval < (0.05 / 85)) {
      print(sprintf("Gene: %s, Outcome: %s: Cauchy p-value = %f", g, p, cauchy_pval))
    }
  }
}



# Function to get Cauchy p-value for gene, phenotype, and tissue
outcome_associated_eqtls_all <- function(gene, outcome) {
  for(tissue in tissues) {
    print(tissue)
  exposure_data <- subset(eqtls, geneSymbol == gene & tissueSiteDetailId == tissue)
  if (nrow(exposure_data) == 0) next
  num_eqtls <- nrow(exposure_data)
  exposure_data <- exposure_data[, c("snpId", "a1", "a2", "beta", "se", "geneSymbol")]
  colnames(exposure_data)[2:5] <- paste0(colnames(exposure_data)[2:5], "_exposure")
  
  outcome_data <- fread(paste0("/Users/kbidding/Documents/mashHF/data/GWAS_subset/", outcome, ".csv"))
  outcome_data$V1 <- NULL
  colnames(outcome_data)[2:5] <- paste0(colnames(outcome_data)[2:5], "_outcome")
  z <- outcome_data$beta_outcome / outcome_data$se_outcome
  outcome_data$pval_outcome <- 2 * pnorm(-abs(z))
  
  merged_data <- merge(exposure_data, outcome_data, by.x = "snpId", by.y = "rsid")
  if (nrow(merged_data) == 0) next
  
  merged_data <- subset(merged_data, (a1_exposure == ea_outcome & a2_exposure == nea_outcome) | 
                          (a2_exposure == ea_outcome & a1_exposure == nea_outcome))
  
  merged_data$beta_outcome[merged_data$a1_exposure == merged_data$nea_outcome] <- 
    -1 * merged_data$beta_outcome[merged_data$a1_exposure == merged_data$nea_outcome]
  
  if (nrow(merged_data) == 0) next
  
  print(merged_data)
  }
}
