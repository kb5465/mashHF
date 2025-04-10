# Load required libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)

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
    labs(title = paste0("Covariance Matrix Sigma_", k), x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
  
  # Print heatmap
  print(p1)
  
  # Save heatmap
  ggsave(paste0("/Users/kbidding/Documents/mashHF/Figures/Sigma_k/heatmap_sigma_", k, ".png"), p1, width = 6, height = 5)
  
  ### Eigenvector barplots ###
  
  # Perform SVD on covariance matrix
  svd_result <- svd(covmat[[k]])
  v <- svd_result$v[h, ]
  d <- svd_result$d
  
  # Loop through first 3 eigenvectors
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
      labs(title = paste("Eigenvector", g, "of Sigma_", k), x = "", y = "Scaled Value") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title.y = element_text(size = 12, color = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "none") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(length(names[h])))
    
    # Save barplot
    ggsave(paste0("/Users/kbidding/Documents/mashHF/Figures/Sigma_k_eigenvectors/U_k", k, "_eigenvector_", g, ".png"), p2, width = 6, height = 4)
    
    # Print barplot
    print(p2)
  }
}