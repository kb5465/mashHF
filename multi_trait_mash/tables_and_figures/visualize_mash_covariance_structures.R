#!/usr/bin/env Rscript

###############################################################
# mashHF: Visualize Covariance Matrices and Eigenvectors
# Author: Kiran Biddinger
# Description:
# - For each U_ed covariance matrix:
#   • Plot heatmap of trait-trait correlation (Σₖ)
#   • Plot first 3 eigenvectors as scaled barplots
#   • Save plots per k and as multi-panel figures
###############################################################

# ----------------------------
# Load Required Libraries
# ----------------------------
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(patchwork)

# ----------------------------
# Initialize Storage Lists
# ----------------------------
heatmap_list <- list()
eigenvector_plots <- list()

# Define trait names and factor levels
trait_order <- c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                 "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                 "NTproBNP", "Troponin")

# ----------------------------
# Loop Over U_ed Matrices
# ----------------------------
for (k in 1:6) {
  
  # Compute correlation matrix
  x <- cov2cor(U_ed[[k]])
  colnames(x) <- trait_order
  rownames(x) <- trait_order
  
  # Melt for ggplot
  lat_melt <- melt(x[rev(trait_order), rev(trait_order)], na.rm = TRUE)
  lat_melt$Var1 <- factor(lat_melt$Var1, levels = rev(trait_order))
  lat_melt$Var2 <- factor(lat_melt$Var2, levels = trait_order)
  
  # Heatmap
  p1 <- ggplot(lat_melt, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", mid = "lightgrey", high = "blue", 
                         midpoint = 0, name = "Correlation", na.value = "grey") +
    labs(title = paste0("Sigma_", k), x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
  
  heatmap_list[[k]] <- p1
  
  # ----------------------------
  # Eigenvector Barplots
  # ----------------------------
  svd_result <- svd(U_ed[[k]])
  v <- svd_result$v
  
  eigen_plots <- list()
  for (g in 1:3) {
    df_eigen <- data.frame(
      Trait = trait_order,
      Eigenvalue = v[, g] / v[which.max(abs(v[, g])), g]
    )
    
    df_eigen$Trait <- factor(df_eigen$Trait, levels = trait_order)
    
    p2 <- ggplot(df_eigen, aes(x = Trait, y = Eigenvalue, fill = Trait)) +
      geom_hline(yintercept = 0, color = "grey") +
      geom_bar(stat = "identity", color = "black", width = 0.9) +
      labs(title = paste("Eigenvector", g), x = "", y = "Scaled Value") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title.y = element_text(size = 12),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "none") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(length(trait_order)))
    
    eigen_plots[[g]] <- p2
  }
  
  # Combine eigenvector plots into single row
  eigenvector_combined <- eigen_plots[[1]] + eigen_plots[[2]] + eigen_plots[[3]] + 
    plot_layout(nrow = 1)
  
  # Save individual eigenvector plot
  ggsave(paste0("/Users/kbidding/Documents/mashHF/Figures/Sigma_k_eigenvectors/U_k", k, "_eigenvectors.png"), 
         eigenvector_combined, width = 12, height = 4)
  
  eigenvector_plots[[k]] <- eigenvector_combined
}

# ----------------------------
# Combine All Heatmaps into a 2x3 Grid
# ----------------------------
heatmap_combined <- heatmap_list[[1]] + heatmap_list[[2]] + heatmap_list[[3]] +
  heatmap_list[[4]] + heatmap_list[[5]] + heatmap_list[[6]] +
  plot_layout(nrow = 2)

# Save the heatmap panel
ggsave("/Users/kbidding/Documents/mashHF/Figures/Sigma_k/heatmap_combined.png", 
       heatmap_combined, width = 12, height = 8)

# Display the combined heatmaps
print(heatmap_combined)
