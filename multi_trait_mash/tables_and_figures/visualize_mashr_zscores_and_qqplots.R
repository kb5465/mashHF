#!/usr/bin/env Rscript

###############################################
# mashHF: Z-score Distribution and QQ Plots
# Author: Kiran Biddinger
# Description:
# - Visualize Z-score density plots for all traits
# - Annotate plots with t-test p-value and AD test p-value
# - Generate QQ plots of p-values for mash_prior and mash_post
###############################################

# ----------------------------
# 1. Load Required Libraries
# ----------------------------
library(ggplot2)
library(cowplot)
library(nortest)

# ----------------------------
# 2. Identify Z-score Columns
# ----------------------------
z_cols <- grep("^Z_", names(mash_prior), value = TRUE)

# ----------------------------
# 3. Create Density Plots
# ----------------------------
plot_list <- list()

for (col_name in z_cols) {
  z_vals <- na.omit(mash_prior[[col_name]])
  
  t_pval <- t.test(z_vals, mu = 0)$p.value
  ad_pval <- ad.test(z_vals)$p.value
  
  p <- ggplot(data.frame(Z = z_vals), aes(x = Z)) +
    geom_density(fill = "steelblue", alpha = 0.3, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.8) +
    labs(
      title = gsub("Z_", "", col_name),
      x = "Z-score",
      y = "Density"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    annotate("text",
             x = Inf, y = Inf,
             hjust = 1.1, vjust = 1.5,
             size = 5,
             label = paste0(
               "t-test p = ", formatC(t_pval, format = "e", digits = 2),
               "\nAD test p = ", formatC(ad_pval, format = "e", digits = 2)
             )
    )
  
  plot_list[[col_name]] <- p
}

# ----------------------------
# 4. Arrange and Save Density Plots
# ----------------------------
top_grid <- plot_grid(plotlist = plot_list[1:12], ncol = 4)
bottom_centered <- ggdraw() +
  draw_plot(plot_list[[13]], x = 0.5, y = 0.5, width = 0.25, height = 1, hjust = 0.5, vjust = 0.5)

final_plot <- plot_grid(top_grid, bottom_centered, ncol = 1, rel_heights = c(3, 1))
print(final_plot)

ggsave("/Users/kbidding/Documents/mashHF/Figures/Zscore_density_plots.pdf", final_plot, width = 14, height = 12)

# ----------------------------
# 5. Define QQ Plot Function
# ----------------------------
plot_qq_pvals <- function(pvals, trait_name) {
  pvals <- na.omit(pvals)
  obs <- sort(pvals)
  exp <- (seq_along(obs)) / (length(obs) + 1)
  
  df <- data.frame(
    observed = -log10(obs),
    expected = -log10(exp)
  )
  
  ggplot(df, aes(x = expected, y = observed)) +
    geom_point(size = 1.4, color = "darkblue", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = trait_name,
      x = expression("Expected " * -log[10](italic(p))),
      y = expression("Observed " * -log[10](italic(p)))
    ) +
    theme_minimal(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# ----------------------------
# 6. Wrapper Function to Create QQ Plots from Dataset
# ----------------------------
create_qqplot_list <- function(df) {
  z_cols <- grep("^Z_", names(df), value = TRUE)
  plot_list <- list()
  
  for (col_name in z_cols) {
    z_vals <- df[[col_name]]
    pvals <- 2 * pnorm(-abs(z_vals))  # two-sided
    short_label <- gsub("^Z_", "", col_name)
    plot_list[[col_name]] <- plot_qq_pvals(pvals, trait_name = short_label)
  }
  
  final_qq <- plot_grid(plotlist = plot_list, ncol = 4)
  return(final_qq)
}

# ----------------------------
# 7. Generate and Save QQ Plots
# ----------------------------
qq_prior <- create_qqplot_list(mash_prior)
print(qq_prior)
ggsave("/Users/kbidding/Documents/mashHF/Figures/qqplots_mash_prior.pdf", qq_prior, width = 14, height = 10)

qq_post <- create_qqplot_list(mash_post)
print(qq_post)
ggsave("/Users/kbidding/Documents/mashHF/Figures/qqplots_mash_post.pdf", qq_post, width = 14, height = 10)
