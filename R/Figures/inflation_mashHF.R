PHENO <- c()
LAMBDA_50 <- c()
LAMBDA_95 <- c()

for (pheno in gsub("BETA_", "", colnames(Bhat))) {
  # Copy dataframe
  dat <- as.data.frame(mash_prior)
  PHENO <- c(PHENO, pheno)
  
  # Calculate p-value from Z-score
  dat$PVAL <- 2 * (1 - pnorm(abs(dat[[paste0("Z_", pheno)]])))

  # Calculate inflation statistic
  chisq <- qchisq(dat$PVAL, lower.tail = FALSE, df = 1)
  
  # Calculate inflation statistic @ 50th percentile (median)
  lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, lower.tail = FALSE, df = 1)
  LAMBDA_50 <- c(LAMBDA_50, lambda_50)
  
  # Write inflation statistic @ 95th percentile
  lambda_95 <- quantile(chisq, 0.95, na.rm = TRUE) / qchisq(0.05, df = 1, lower.tail = FALSE)
  LAMBDA_95 <- c(LAMBDA_95, lambda_95)
  
  # Print results
  print(paste0(pheno, ": lambda_50 = ", lambda_50, "; lambda_95 = ", lambda_95))
}

# Make data frame of results
results <- data.frame(PHENO, LAMBDA_50, LAMBDA_95)
