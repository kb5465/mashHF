# Load required libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(patchwork) # For multi-panel plots

heatmap_list <- list()  # Store individual heatmaps
eigenvector_plots <- list()  # Store eigenvector plots per k

# Loop over k from 1 to 6
for (k in 1:6) {
  
  # Compute the correlation matrix
  covmat <- U_ed
  x <- cov2cor(covmat[[k]])
  names <- gsub("BETA_", "", colnames(Bhat)) # Trait names
  colnames(x) <- names
  rownames(x) <- names
  
  # Indices of outcome traits (all in this case)
  h <- rep(TRUE, length(names))
  names <- factor(names, levels = c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                                    "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                                    "NTproBNP", "Troponin"))
  
  # Convert matrix to long format for ggplot
  lat <- x[rev(h), rev(h)]
  lat_melt <- melt(lat, na.rm = TRUE)
  lat_melt$Var1 <- factor(lat_melt$Var1, levels = rev(c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                                                        "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                                                        "NTproBNP", "Troponin")))
  lat_melt$Var2 <- factor(lat_melt$Var2, levels = c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                                                    "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                                                    "NTproBNP", "Troponin"))
  
  # Set color palette for heatmap
  mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(100)
  
  # Generate heatmap using ggplot2
  p1 <- ggplot(lat_melt, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", mid = "lightgrey", high = "blue", midpoint = 0, name = "Correlation", na.value = "grey") +
    labs(title = paste0("Sigma_", k), x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
  
  # Store heatmap
  heatmap_list[[k]] <- p1
  
  ### Eigenvector barplots ###
  
  # Perform SVD on covariance matrix
  svd_result <- svd(covmat[[k]])
  v <- svd_result$v[h, ]
  d <- svd_result$d
  
  # Store eigenvector plots
  eigen_plots <- list()
  
  for (g in 1:3) {
    
    # Create a data frame for ggplot
    df_eigen <- data.frame(
      Trait = names[h],
      Eigenvalue = v[, g] / v[which.max(abs(v[, g])), g]
    )
    
    # Factor levels
    df_eigen$Trait <- factor(df_eigen$Trait, levels = c("HF", "DCM", "HCM", "LVEF", "LVESVi", "LVEDVi",
                                                        "LVMi", "meanWT", "strain_circ", "strain_long", "strain_rad", 
                                                        "NTproBNP", "Troponin"))
    
    # Generate barplot with Blues color scheme
    p2 <- ggplot(df_eigen, aes(x = Trait, y = Eigenvalue, fill = Trait)) +
      geom_hline(yintercept = 0, color = "grey") +
      geom_bar(stat = "identity", color = "black", width = 0.9) +
      labs(title = paste("Eigenvector", g), x = "", y = "Scaled Value") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title.y = element_text(size = 12, color = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "none") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(length(names[h])))
    
    # Store eigenvector plot
    eigen_plots[[g]] <- p2
  }
  
  # Combine three eigenvector plots into one row
  eigenvector_combined <- eigen_plots[[1]] + eigen_plots[[2]] + eigen_plots[[3]] + 
    plot_layout(nrow = 1)
  
  # Save the combined eigenvector plot
  ggsave(paste0("/Users/kbidding/Documents/mashHF/Figures/Sigma_k_eigenvectors/U_k", k, "_eigenvectors.png"), 
         eigenvector_combined, width = 12, height = 4)
  
  # Store the combined eigenvector plot
  eigenvector_plots[[k]] <- eigenvector_combined
}

# Combine all heatmaps into one figure (3 per row)
heatmap_combined <- heatmap_list[[1]] + heatmap_list[[2]] + heatmap_list[[3]] +
  heatmap_list[[4]] + heatmap_list[[5]] + heatmap_list[[6]] + 
  plot_layout(nrow = 2)

# Save the combined heatmap plot
ggsave("/Users/kbidding/Documents/mashHF/Figures/Sigma_k/heatmap_combined.png", heatmap_combined, width = 12, height = 8)

# Print the combined heatmap plot
print(heatmap_combined)
