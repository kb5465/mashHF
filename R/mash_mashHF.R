library(data.table)
library(mashr)
library(tidyverse)
library(biomaRt)
library(reshape2)
library(gaston)

# Read in merged Beta and SE statistics
rvat <- fread("/Users/kbidding/Documents/mashHF/data/mashHF_canonical_LOF_0.01_merged_firth.tsv")

# Remove genes with NA tests
mash_prior <- drop_na(rvat)
mash_prior <- as.data.frame(mash_prior)

# Subset to phenotypes
# outcomes <- c("DCM", "HCM")
# columns_to_keep <- c("ID", paste0("BETA_", outcomes), paste0("SE_", outcomes))
# mash_prior <- mash_prior[, columns_to_keep]

# Extract Bhat and Shat matrices
Bhat <- as.matrix(mash_prior[, grep("^BETA_", colnames(mash_prior))])
Shat <- as.matrix(mash_prior[, grep("^SE_", colnames(mash_prior))])

# Prepare the data object
mash_data <- mash_set_data(Bhat, Shat)

# Calculate z-statistics for each phenotype
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  mash_prior[[paste0("Z_", pheno)]] <- mash_prior[[paste0("BETA_", pheno)]] / mash_prior[[paste0("SE_", pheno)]]
}

# Identify a subset of strong tests
m_1by1 <- mash_1by1(mash_data)

# Estimate null correlation structure
Vhat <- estimate_null_correlation_simple(mash_data)
image(Vhat)

strong_subset <- get_significant_results(m_1by1, 0.1) ## N=49 strong genes
set.seed(123)
remaining_indices <- setdiff(1:nrow(Bhat), strong_subset)
random_subset <- sample(remaining_indices, 100)

# Extract subset of strong tests
mash_strong <- mash_set_data(Bhat[strong_subset,], Shat[strong_subset,], V=Vhat)
mash_random <- mash_set_data(Bhat[random_subset,], Shat[random_subset,], V=Vhat)
# mash_strong <- mash_set_data(Bhat[strong_subset,], Shat[strong_subset,], V=diag(nrow = nrow(Vhat)))

# Create data-driven covariance matrices
U_pca <- cov_pca(mash_strong, 5)
U_ed <- cov_ed(mash_strong, U_pca)
# Create canonical covariance matrices
U_cov <- cov_canonical(mash_data)

g=1
# Identify the max absolute value in column g of U_ed
max_value <- max(abs(U_ed[[g]]), na.rm = TRUE)

# Normalize the eigenvector by its maximum absolute value
normalized_vector <- U_ed[[g]] / max_value

# Create the bar plot
barplot(normalized_vector, 
        las = 2, 
        main = paste("Eigenvector", g, "of Uk"), 
        cex.names = 0.5)

# Run mash
mash_result <- mash(mash_data, Ulist = c(U_ed, U_cov), pi_thresh = 1e-30)

qqplot.pvalues(as.data.frame(mash_result$result$lfsr)$BETA_DCM)

plot(density(mash_prior$Z_HCM))
plot(density(as.data.frame(mash_result$result$PosteriorMean)$BETA_DCM))

# Print likelihood weights
print(get_estimated_pi(mash_result))

# Get inflation statistics
get_lambda_singletrait(B_obs, SE)
get_lambda_post(mash_result)

# Extract posterior means
mash_post_mean <- as.data.frame(mash_result$result$PosteriorMean)
mash_post_sd <- as.data.frame(mash_result$result$PosteriorSD)
mash_post_lfsr <- as.data.frame(mash_result$result$lfsr)

# Rename rows to ENSEMBL gene IDs
mash_post_mean$ID <- mash_prior$ID
mash_post_sd$ID <- mash_prior$ID
mash_post_lfsr$ID <- mash_prior$ID
colnames(mash_post_sd) <- gsub("BETA", "SE", colnames(mash_post_sd))
colnames(mash_post_lfsr) <- gsub("BETA", "lfsr", colnames(mash_post_lfsr))

# Make a mash posterior data frame
mash_post <- as.data.frame(merge(mash_post_mean, mash_post_sd, by = "ID"))
mash_post <- merge(mash_post, mash_post_lfsr, by = "ID")

# Calculate z-statistics for each phenotype
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  mash_prior[[paste0("Z_", pheno)]] <- mash_prior[[paste0("BETA_", pheno)]] / mash_prior[[paste0("SE_", pheno)]]
  mash_post[[paste0("Z_", pheno)]] <- mash_post[[paste0("BETA_", pheno)]] / mash_post[[paste0("SE_", pheno)]]
}

# Compare previous and post
subset(rvat, rvat$ID == "ENSG00000155657")
subset(mash_post_mean, mash_post_mean$ID == "ENSG00000155657")
subset(mash_post_sd, mash_post_sd$ID == "ENSG00000155657")

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch gene annotations
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = mash_post$ID,
  mart = ensembl
)

# Rename columns for merging
colnames(gene_annotations) <- c("ENSEMBL_ID", "Gene_Name")

# Merge with mash_post
mash_post_annotated <- mash_post %>%
  left_join(gene_annotations, by = c("ID" = "ENSEMBL_ID"))

# Merge with mash_prior
mash_prior_annotated <- mash_prior %>%
  left_join(gene_annotations, by = c("ID" = "ENSEMBL_ID"))

# View the merged data
head(mash_post_annotated)


# Order by most significant results for a given phenotype
mash_post_annotated %>%
  arrange(desc(abs(Z_HF))) %>% 
  slice(1:10) %>% 
  pull(Gene_Name)

# Calculate the Bonferroni significance threshold
bonf_sig_thr <- 0.05 / nrow(mash_prior)

# Calculate the corresponding z-statistic threshold for two-tailed test
z_thr <- qnorm(1 - bonf_sig_thr / 2)

# Pull exome-wide significant hits per phenotype -- PRIOR
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  # Subset to exome-wide significant hits
  sig_hits <- subset(mash_prior_annotated, abs(mash_prior_annotated[[paste0("Z_", pheno)]]) > z_thr)
  
  # Calculate p-value from standard errors
  sig_hits$PVAL <- 2 * (1 - pnorm(abs(sig_hits[[paste0("Z_", pheno)]])))
  
  # Convert to data frame
  sig_hits <- as.data.frame(sig_hits)
  
  # Subset to columns of interest
  sig_hits <- sig_hits[,c("Gene_Name", paste0("BETA_", pheno), paste0("SE_", pheno), "PVAL")]
  
  print(sig_hits)
  
  # Write tables
  # write.csv(sig_hits, paste0("/Users/kbidding/Documents/mashHF/data/", pheno, "_prior_sighits.csv"))
}

# Pull exome-wide significant hits per phenotype -- POSTERIOR
bonf_sig_thr <- 0.05 / nrow(mash_post)
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  # Subset to exome-wide significant hits
  sig_hits <- subset(mash_post_annotated, abs(mash_post_annotated[[paste0("Z_", pheno)]]) > z_thr)
  
  # Calculate p-value from Z-score
  sig_hits$PVAL <- 2 * (1 - pnorm(abs(sig_hits[[paste0("Z_", pheno)]])))
  
  # Subset to columns of interest
  sig_hits <- sig_hits[,c("Gene_Name", paste0("BETA_", pheno), paste0("SE_", pheno), "PVAL", paste0("lfsr_", pheno))]
  
  # Write tables 
  print(pheno)
  print(dim(sig_hits))
  print("CHI3L1" %in% sig_hits$Gene_Name)
  print(subset(sig_hits, Gene_Name == "CHI3L1"))
  print(head(sig_hits))
}

# Compare pre and post
# Define the highlighted genes
genes_highlight <- c("TTN", "BAG3", "DES", "FLNC", "LMNA", "MYH7", "PLN", "RBM20", 
                     "SCN5A", "TNNC1", "TNNT2", "DSP",
                     "MYBPC3")
# Iterate over each gene
for (gene in genes_highlight) {
  # Subset prior and posterior datasets for the given gene
  prior <- as.data.frame(subset(mash_prior_annotated, Gene_Name == gene))
  posterior <- as.data.frame(subset(mash_post_annotated, Gene_Name == gene))
  
  # If empty data, skip
  if (nrow(prior) == 0 | nrow(posterior) == 0) {
    cat("\nSkipping gene:", gene, "due to no prior or posterior data.\n")
    next
  }
  
  # Select only relevant column names that contain BETA, SE, or Z
  prior_selected <- names(prior)[grepl("BETA|SE|Z", names(prior))]
  posterior_selected <- names(posterior)[grepl("BETA|SE|Z", names(posterior))]
  
  # Find common metrics between prior and posterior datasets
  common_metrics <- intersect(prior_selected, posterior_selected)
  
  # If no common metrics, skip the gene
  if (length(common_metrics) == 0) {
    cat("\nSkipping gene:", gene, "due to no common metrics.\n")
    next
  }
  
  # Subset only common metrics
  prior_selected <- prior[, common_metrics]
  posterior_selected <- posterior[, common_metrics]
  
  # Combine prior and posterior by row and transpose
  combined_df <- rbind(prior_selected, posterior_selected)
  transposed_df <- as.data.frame(t(combined_df))
  
  # Rename columns for clarity
  colnames(transposed_df) <- c("Prior", "Posterior")
  transposed_df$Metric <- rownames(transposed_df)
  
  # Print the dataframe for the current gene
  cat("\nComparison for gene:", gene, "\n")
  print(transposed_df)
}



#### HELPER FUNCTIONS #########
# Helper function to get lambda from observed single-trait dada
get_lambda_singletrait <- function(B_obs, SE) {
  n_col <- ncol(B_obs)
  results <- data.frame()
  
  # Iterate over columns
  for (i in 1:n_col) {
    # Pull columns for specific traits
    beta = as.data.frame(B_obs)[, i]
    se = as.data.frame(SE)[, i]
    
    # Get z-scores
    z <- beta / se
    
    # Extract p-values 
    pval <- 2 * pnorm(-abs(z), mean = 0, sd = 1)
    
    # Compute chi-square statistics
    chisq <- z^2
    
    # OPTIONAL: Print QQ plot
    print(qqplot.pvalues(pval, main = paste0("ST: Trait_", i)))
    print(density(pval))
    
    # Compute genomic inflation factor at 50th percentile
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    
    # Add data to list
    results <- rbind(results, data.frame(Trait = paste0("Trait_", i), Inflation = lambda_50))
  }
  
  # Get z-score
  # Add column names
  colnames(results) = c("Trait", "Inflation")
  results$Inflation <- as.numeric(results$Inflation)
  
  return(results)
}

# Helper function to get lambda from mash result
get_lambda_post <- function(data) {
  n_col <- ncol(data$result$PosteriorMean)
  
  colnames(data$result$PosteriorMean) <- paste0("BETA_", seq(1, n_col))
  colnames(data$result$PosteriorSD) <- paste0("SE_", seq(1, n_col))
  
  # Set to dataframes
  data$result$PosteriorMean <- as.data.frame(data$result$PosteriorMean)
  data$result$PosteriorSD <- as.data.frame(data$result$PosteriorSD)
  
  results <- data.frame()
  
  # Iterate over columns
  for (i in 1:n_col) {
    # Pull columns for specific traits
    beta = data$result$PosteriorMean[[paste0("BETA_", i)]]
    se = data$result$PosteriorSD[[paste0("SE_", i)]]
    
    # Get z-scores
    z <- beta / se
    
    # Extract p-values 
    pval <- 2 * pnorm(-abs(z), mean = 0, sd = 1)
    
    # Compute chi-square statistics
    chisq <- z^2
    
    # OPTIONAL: Print QQ plot
    print(qqplot.pvalues(pval, main = paste0("MT: Trait_", i)))
    
    # Compute genomic inflation factor at 50th percentile
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    
    # Add data to list
    results <- rbind(results, data.frame(Trait = paste0("Trait_", i), Inflation = lambda_50))
  }
  
  # Get z-score
  # Add column names
  colnames(results) = c("Trait", "Inflation")
  results$Inflation <- as.numeric(results$Inflation)
  
  return(results)
}




### ---- WRITE TABLES --- ####
library(openxlsx)

# Define Bonferroni threshold
bonf_sig_thr <- 0.05 / nrow(mash_prior)
z_thr <- qnorm(1 - bonf_sig_thr / 2)

# Initialize an Excel workbook
wb <- createWorkbook()

# Iterate through phenotypes
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  
  # Subset to exome-wide significant hits
  sig_hits <- subset(mash_prior_annotated, abs(mash_prior_annotated[[paste0("Z_", pheno)]]) > z_thr)
  
  # Calculate p-value from Z-score
  sig_hits$PVAL <- 2 * (1 - pnorm(abs(sig_hits[[paste0("Z_", pheno)]])))
  
  # Make data frame
  sig_hits <- as.data.frame(sig_hits)
  
  # Subset to columns of interest
  sig_hits <- sig_hits[, c("Gene_Name", paste0("BETA_", pheno), paste0("SE_", pheno), "PVAL")]
  
  # Add sheet and write data
  addWorksheet(wb, pheno)
  writeData(wb, pheno, sig_hits)
  
  # Print debug info
  print(pheno)
  print(dim(sig_hits))
  print(head(sig_hits))
}

# Save the Excel file
saveWorkbook(wb, "/Users/kbidding/Documents/mashHF/Tables/mash_prior_significant_hits.xlsx", overwrite = TRUE)





bonf_sig_thr <- 0.05 / nrow(mash_prior)
z_thr <- qnorm(1 - bonf_sig_thr / 2)

# Initialize an Excel workbook
wb <- createWorkbook()

# Iterate through phenotypes
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  
  # Subset to exome-wide significant hits
  sig_hits <- subset(mash_post_annotated, abs(mash_post_annotated[[paste0("Z_", pheno)]]) > z_thr)
  
  # Calculate p-value from Z-score
  sig_hits$PVAL <- 2 * (1 - pnorm(abs(sig_hits[[paste0("Z_", pheno)]])))
  
  # Subset to columns of interest
  sig_hits <- sig_hits[, c("Gene_Name", paste0("BETA_", pheno), paste0("SE_", pheno), "PVAL")]
  
  # Add sheet and write data
  addWorksheet(wb, pheno)
  writeData(wb, pheno, sig_hits)
  
  # Print debug info
  print(pheno)
  print(dim(sig_hits))
  print(head(sig_hits))
}

# Save the Excel file
saveWorkbook(wb, "/Users/kbidding/Documents/mashHF/Tables/mash_post_significant_hits.xlsx", overwrite = TRUE)





