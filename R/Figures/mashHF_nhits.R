library(ggplot2)
library(dplyr)
library(data.table)
library(tibble)
library(tidyr)
library(RColorBrewer)

# Helper function to get lambda from data
get_nhits <- function(data) {
  results <- data.frame()
  
  # Iterate over traits to create outcome dataframe
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    # Extract p-values 
    pval <- 2 * pnorm(abs(data[[paste0("Z_", pheno)]]), mean = 0, sd = 1, lower.tail = F)
    
    # Calculate Bonferonni significance level
    bonf <- 0.05 / sum(!is.na(pval))

    # Calculate the number of significant hits
    prop_hits <- sum(pval < bonf) #/ sum(!is.na(pval)) * 100
    
    # Add data to list
    results <- rbind(results, c(pheno, prop_hits))
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
barplot_nhits_fn <- function(data, plotname, title) {
  # Get results using data
  results <- get_nhits(data)
  
  # Set color palette
  mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(13)
  
  # Barplot of results
  p <- ggplot(results, aes(x = Trait, y = Inflation)) +
    geom_hline(yintercept = 0, color = "grey") + 
    geom_bar(aes(fill = Trait), stat = "identity", color = "black", width = 0.9) +
    labs(x = "", y = "Number of Genes Associated", title = title) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1), 
          axis.title.x = element_text(size = 11, color = "black"),  
          axis.text.y = element_text(size = 11, color = "black"), 
          axis.title.y = element_text(size = 11, color = "black"),
          axis.ticks.y = element_line(color = "black"),
          axis.ticks.x = element_line(color = "black"),
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_fill_manual(values = mycolors)

  # Save plot
  ggsave(plotname, p)
  
  # Print results
  print(p)
}

# Write barplots
barplot_nhits_fn(mash_prior_annotated, "/Users/kbidding/Documents/mashHF/Figures/mash_prior_nhits.png", "Single-Trait")
barplot_nhits_fn(mash_post_annotated, "/Users/kbidding/Documents/mashHF/Figures/mash_post_nhits.png", "Multi-Trait")


