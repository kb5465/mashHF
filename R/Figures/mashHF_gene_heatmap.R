library(ggplot2)
library(dplyr)
library(data.table)
library(tibble)
library(tidyr)
library(RColorBrewer)

# Extract genes associated with >=1 trait
get_all_genes <- function(data) {
  results <- c()
  
  # Iterate over traits to create outcome dataframe
  for (pheno in gsub("BETA_", "", colnames(Bhat))) {
    # Extract p-values 
    pval <- 2 * pnorm(abs(data[[paste0("Z_", pheno)]]), mean = 0, sd = 1, lower.tail = F)
    
    # Calculate Bonferonni significance level
    bonf <- 0.05 / sum(!is.na(pval))
    
    # Extract significant hits
    genes <- data$Gene_Name[pval < bonf] 
    
    # Add data to list
    results <- c(results, genes)
  }
  
  # Subset to unique list
  results <- unique(results)
  results <- results[results != ""]
  
  return(results)
}

# Extract associations with all traits
get_gene_associations <- function(data) {
  # Get significant genes
  genes <- get_all_genes(data)
  
  # Subset data to relevant genes
  data <- subset(data, Gene_Name %in% genes)
  
  # Subset to gene name and z-score columns
  cols <- (colnames(data) == "Gene_Name") | (grepl("Z_", colnames(data)))
  data <- data[, cols]
  
  # Rename columns to get rid of "Z_" prefix
  colnames(data) <- gsub("Z_", "", colnames(data))
  
  # P
  return(data)
}



# Rename columns for clarity
colnames(heatmap_data) <- c("Gene", "Trait", "Association")

# Create the heatmap
ggplot(heatmap_data, aes(x = Trait, y = Gene, fill = Association)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(x = "Outcome Trait", y = "Gene", fill = "Association Strength") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


# Heatmap helper function
genes_heatmap_fn <- function(data, plotname) {
  # Convert data from wide to long format for ggplot
  dat <- get_gene_associations(data)

  # Convert data from wide to long format
  heatmap_data <- dat %>%
    pivot_longer(cols = -Gene_Name, names_to = "Trait", values_to = "Association")
  
  # Factor levels
  heatmap_data$Gene_Name <- factor(heatmap_data$Gene_Name, levels = rev(sort(unique(heatmap_data$Gene_Name))))
  heatmap_data$Trait <- factor(heatmap_data$Trait, levels = c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                                                    "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                                                    "NTproBNP", "Troponin"))
  
  # Create the heatmap
  p <- ggplot(heatmap_data, aes(x = Trait, y = Gene_Name, fill = Association)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    labs(x = "", y = "", fill = "Association Strength (Z)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"))
  
  # Save plot
  ggsave(plotname, p, width = 6.56, height = 12)
  
  # Print results
  print(p)
}

genes_heatmap_fn(mash_post_annotated, "/Users/kbidding/Documents/mashHF/Figures/mash_post_heatmap.png")







library(ggplot2)
library(tidyr)
library(dplyr)

# Convert data from wide to long format
dat <- get_gene_associations(mash_post_annotated)

# Convert to long format
heatmap_data <- dat %>%
  pivot_longer(cols = -Gene_Name, names_to = "Trait", values_to = "Association")

# Sort Gene_Name alphabetically from top to bottom
heatmap_data$Gene_Name <- factor(heatmap_data$Gene_Name, levels = rev(sort(unique(heatmap_data$Gene_Name))))

# Define trait groups
left_traits <- c("HF", "DCM", "HCM")
right_traits <- c("LVEF", "LVESVi", "LVEDVi", "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", "NTproBNP", "Troponin")

# Assign a "Side" to each trait to separate them
heatmap_data$Side <- ifelse(heatmap_data$Trait %in% left_traits, "Left", "Right")

# Order Trait factor so left-side traits appear first
heatmap_data$Trait <- factor(heatmap_data$Trait, levels = c(left_traits, right_traits))

# Create the heatmap
ggplot(heatmap_data, aes(x = Trait, y = Gene_Name, fill = Association)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(x = "", y = "", fill = "Association Strength (Z)") +
  facet_grid(. ~ Side, space = "free", scales = "free", switch = "y", axis.labels = "all_y") +  # Split into two panels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.title.y = element_blank(),  # Remove y-axis title
        strip.text.y.left = element_text(size = 10, angle = 90),  # Move gene names between panels
        strip.placement = "outside",  # Ensure proper placement of labels
        panel.spacing.x = unit(0.5, "cm"))   # Add a small space between left & right panels


