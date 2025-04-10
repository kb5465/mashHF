#===============================#
# Load libraries
#===============================#
library(ggplot2)
library(cowplot)
library(nortest)

#===============================#
# Identify columns of interest
#===============================#
# We'll assume your data frame is called mash_prior.
# We detect all columns whose names start with "Z_":
z_cols <- grep("^Z_", names(mash_prior), value = TRUE)

#===============================#
# Create individual density plots
#===============================#
plot_list <- list()

for (col_name in z_cols) {
  # Extract the vector of Z-scores
  z_vals <- mash_prior[[col_name]]
  
  # Remove NA entries (if any)
  z_vals <- z_vals[!is.na(z_vals)]
  
  # -- 1) t-test for mean != 0
  t_res <- t.test(z_vals, mu = 0)
  t_pval <- t_res$p.value
  
  # -- 2) Test of normality (Anderson-Darling)
  ad_res <- ad.test(z_vals)
  ad_pval <- ad_res$p.value
  
  # -- 3) Build the density plot
  p <- ggplot(data.frame(Z = z_vals), aes(x = Z)) +
    geom_density(fill = "steelblue", alpha = 0.3, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.8) +
    labs(
      title = gsub("Z_", "", col_name),
      x = "Z-score",
      y = "Density"
    ) +
    # Use a larger base font size for figure-readiness
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    # Annotate p-values in top-right corner
    annotate(
      "text", 
      x = Inf, 
      y = Inf, 
      hjust = 1.1,        # shift text slightly left of the right border
      vjust = 1.5,        # shift text slightly below the top border
      size = 5,
      label = paste0(
        "t-test p = ", formatC(t_pval, format = "e", digits = 2),
        "\nAD test p = ", formatC(ad_pval, format = "e", digits = 2)
      )
    )
  
  plot_list[[col_name]] <- p
}

#===============================#
# Arrange all density plots
#===============================#
# 1) Arrange the first 12 plots in a 3-row × 4-column grid
top_grid <- plot_grid(
  plotlist = plot_list[1:12],
  ncol = 4
)

# 2) Center the 13th plot in its own row
#    We'll use `ggdraw()` + `draw_plot()` to precisely place the single plot in the center.
bottom_centered <- ggdraw() + 
  draw_plot(plot_list[[13]], 
            x = 0.5,      # place center of the plot at 50% across the width
            y = 0.5,      # place center of the plot at 50% down the height
            width = 0.25,  # make it as wide as the entire region
            height = 1, # and as tall as the entire region
            hjust = 0.5,  # interpret x as the center, not left edge
            vjust = 0.5)  # interpret y as the center, not bottom edge

# 3) Stack the top grid + the centered bottom plot into one figure
final_plot <- plot_grid(
  top_grid,
  bottom_centered,
  ncol = 1,           
  rel_heights = c(3, 1)  # adjust as needed so the bottom row is shorter
)

# 4) Display or save
print(final_plot)

# Save to file:
ggsave("/Users/kbidding/Documents/mashHF/Figures/Zscore_density_plots.pdf", final_plot, width = 14, height = 12)


#=================================================
# Define a Function to Build a P-value QQ-plot
#=================================================
plot_qq_pvals <- function(pvals, trait_name) {
  # Remove NA
  pvals <- pvals[!is.na(pvals)]
  
  # Sort observed p-values
  obs <- sort(pvals)
  
  # Theoretical uniform quantiles
  exp <- (seq_along(obs)) / (length(obs) + 1)
  
  # Build data frame for ggplot
  df <- data.frame(
    observed = -log10(obs),
    expected = -log10(exp)
  )
  
  # Generate plot
  p <- ggplot(df, aes(x = expected, y = observed)) +
    geom_point(size = 1.4, color = "darkblue", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = trait_name, 
      x = expression("Expected " * -log[10](italic(p))),
      y = expression("Observed " * -log[10](italic(p)))
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

#=================================================
# Function to Create a QQ-plot List from a Dataset
#=================================================
# This expects that `df` has columns like Z_DCM, Z_HCM, etc.
create_qqplot_list <- function(df, dataset_label = "prior") {
  # Identify Z-score columns
  z_cols <- grep("^Z_", names(df), value = TRUE)
  
  # Prepare an empty list to store ggplot objects
  plot_list <- list()
  
  for (col_name in z_cols) {
    # Compute p-values
    z_vals <- df[[col_name]]
    pvals <- 2 * pnorm(-abs(z_vals))  # two-sided
    # Generate a short label by removing "Z_"
    short_label <- gsub("^Z_", "", col_name)
    
    # Create the QQ-plot
    p <- plot_qq_pvals(pvals, trait_name = short_label)
    plot_list[[col_name]] <- p
  }
  
  # Combine them into a grid
  # E.g., use 4 columns to get ~3-4 rows if you have ~13 traits
  final_qq <- plot_grid(plotlist = plot_list, ncol = 4)
  
  return(final_qq)
}

#=================================================
# Create and Plot for mash_prior
#=================================================
qq_prior <- create_qqplot_list(mash_prior)
print(qq_prior)  # to display in your R session

#=================================================
# Create and Plot for mash_post
#=================================================
# If mash_post has columns like Z_DCM, Z_HCM, etc. (the posterior z-scores)
qq_post <- create_qqplot_list(mash_post)
print(qq_post)

#=================================================
# (Optional) Save to Files
#=================================================
ggsave("/Users/kbidding/Documents/mashHF/Figures/qqplots_mash_prior.pdf", qq_prior, width = 14, height = 10)
ggsave("/Users/kbidding/Documents/mashHF/Figures/qqplots_mash_post.pdf", qq_post, width = 14, height = 10)

