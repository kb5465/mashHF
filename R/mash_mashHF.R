library(data.table)
library(mashr)
library(tidyverse)
library(dplyr)

# Read in merged Beta and SE statistics
rvat <- fread("/Users/kbidding/Documents/mashHF/data/mashHF_canonical_LOF_0.01_merged.tsv")

# Remove genes with NA tests
rvat_cleaned <- drop_na(rvat)

# Extract Bhat and Shat matrices
Bhat <- as.matrix(rvat_cleaned[, grep("^BETA_", colnames(rvat)), with = FALSE])
Shat <- as.matrix(rvat_cleaned[, grep("^SE_", colnames(rvat)), with = FALSE])

# Prepare the data object
mash_data <- mash_set_data(Bhat, Shat)

# Identify a subset of strong tests
m_1by1 <- mash_1by1(mash_data)
strong_subset <- get_significant_results(m_1by1, 0.05)

# Estimate null correlation structure
Vhat <- estimate_null_correlation_simple(mash_data)

# Extract subset of strong tests
mash_strong <- mash_set_data(Bhat[strong_subset,], Shat[strong_subset,], V=Vhat)

# Create data-driven covariance matrices
U_pca <- cov_pca(mash_strong,5)
U_ed <- cov_ed(mash_strong, U_pca)

# Create canonical covariance matrices
U_cov <- cov_canonical(mash_data)

# Run mash
mash_result <- mash(mash_data, Ulist = c(U_ed, U_cov))

# Extract posterior means
mash_post_mean <- as.data.frame(mash_result$result$PosteriorMean)
mash_post_sd <- as.data.frame(mash_result$result$PosteriorSD)

# Rename rows to ENSEMBL gene IDs
mash_post_mean$ID <- rvat_cleaned$ID
mash_post_sd$ID <- rvat_cleaned$ID
colnames(mash_post_sd) <- gsub("BETA", "SE", colnames(mash_post_sd))

# Make a mash posterior data frame
mash_post <- as.data.frame(merge(mash_post_mean, mash_post_sd, by = "ID"))

# Calculate z-statistics for each phenotype
for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  mash_post[[paste0("Z_", pheno)]] <- mash_post[[paste0("BETA_", pheno)]] / mash_post[[paste0("SE_", pheno)]]
}

# Order by most significant results for a given phenotype
mash_post %>%
  arrange(desc(abs(Z_DCM))) %>% 
  slice(1:10)


# Compare TTN gene row
subset(rvat, rvat$ID == "ENSG00000155657")
subset(mash_post_mean, mash_post_mean$ID == "ENSG00000155657")
subset(mash_post_sd, mash_post_sd$ID == "ENSG00000155657")

