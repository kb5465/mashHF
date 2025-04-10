library(ggplot2)
library(dplyr)
library(data.table)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(gaston)

# Helper function to get lambda from data
get_lambda <- function(data) {
  # Initialize empty dataframe
  results <- data.frame()
  
  # Iterate over traits to create outcome dataframe
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    # Extract p-values 
    pval <- 2 * pnorm(-abs(data[[paste0("Z_", pheno)]]), mean = 0, sd = 1)
    
    # Compute chi-square statistics
    chisq <- data[[paste0("Z_", pheno)]]^2
    
    # OPTIONAL: Print QQ plot
    print(qqplot.pvalues(pval, main = pheno, ylim = c(0, 10)))
    
    # Compute genomic inflation factor at 50th percentile
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    
    # Add data to list
    results <- rbind(results, data.frame(Trait = pheno, Inflation = lambda_50))
  }
  
  # Add column names
  colnames(results) = c("Trait", "Inflation")
  results$Inflation <- as.numeric(results$Inflation)
  
  # Factor ordering for plot
  results$Trait <- factor(results$Trait, levels = c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                                                    "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                                                    "NTproBNP", "Troponin"))
  
  return(results)
}

# Barplot helper function
barplot_fn <- function(data, plotname, title) {
  # Get results using data
  results <- get_lambda(data)

  # Set color palette
  mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(13)
  
  # Barplot of results
  p <- ggplot(results, aes(x = Trait, y = Inflation)) +
    geom_hline(yintercept = 1, color = "grey") + 
    geom_bar(aes(fill = Trait), stat = "identity", color = "black", width = 0.9) +
    labs(x = "", y = "Genomic Inflation", title = title) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1), 
          axis.title.x = element_text(size = 11, color = "black"),  
          axis.text.y = element_text(size = 11, color = "black"), 
          axis.title.y = element_text(size = 11, color = "black"),
          axis.ticks.y = element_line(color = "black"),
          axis.ticks.x = element_line(color = "black"),
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_fill_manual(values = mycolors) +
    ylim(0, 2.3)
  
  # Save plot
  ggsave(plotname, p)
  
  # Print results
  print(p)
}

# Write barplots
barplot_fn(mash_prior_annotated, "/Users/kbidding/Documents/mashHF/Figures/mash_prior_lambda.png", "Single-Trait")
barplot_fn(mash_post_annotated, "/Users/kbidding/Documents/mashHF/Figures/mash_post_lambda.png", "Multi-Trait")






# Helper function to get lambda at 95th percentile from data
get_lambda_95 <- function(data) {
  results <- data.frame()
  
  # Iterate over traits to create outcome dataframe
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    # Extract p-values 
    pval <- 2 * pnorm(abs(data[[paste0("Z_", pheno)]]), mean = 0, sd = 1, lower.tail = F)
    
    # Calculate inflation statistic
    chisq <- qchisq(pval, lower.tail = FALSE, df = 1)
    
    # Calculate inflation statistic @ 95th percentile
    lambda_95 <- quantile(chisq, 0.95, na.rm = TRUE) / qchisq(0.05, df = 1, lower.tail = FALSE)
    
    # Add data to list
    results <- rbind(results, c(pheno, lambda_95))
  }
  
  # Add column names
  colnames(results) = c("Trait", "Inflation")
  results$Inflation <- as.numeric(results$Inflation)
  
  # Factor ordering for plot
  results$Trait <- factor(results$Trait, levels = c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                                                    "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                                                    "NTproBNP", "Troponin"))
  
  return(results)
}

# Helper function to get lambda at 95th percentile from data
get_lambda_5 <- function(data) {
  results <- data.frame()
  
  # Iterate over traits to create outcome dataframe
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    # Extract p-values 
    pval <- 2 * pnorm(abs(data[[paste0("Z_", pheno)]]), mean = 0, sd = 1, lower.tail = F)
    
    # Calculate inflation statistic
    chisq <- qchisq(pval, lower.tail = FALSE, df = 1)
    
    # Calculate inflation statistic @ 5th percentile
    lambda_5 <- quantile(chisq, 0.05, na.rm = TRUE) / qchisq(0.95, df = 1, lower.tail = FALSE)
    
    # Add data to list
    results <- rbind(results, c(pheno, lambda_5))
  }
  
  # Add column names
  colnames(results) = c("Trait", "Inflation")
  results$Inflation <- as.numeric(results$Inflation)
  
  # Factor ordering for plot
  results$Trait <- factor(results$Trait, levels = c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                                                    "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                                                    "NTproBNP", "Troponin"))
  
  return(results)
}