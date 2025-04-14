########################################################################
# Simulation: Compare Effect Structures Under Fixed Settings
# Settings:
# - Structured noise
# - 1% genes with real effects
# - Weak trait correlation (cov1 = 0.3, cov2 = 0.7)
# Varying:
# - Shared vs. Specific vs. Patterned effects
########################################################################

# Define effect types to iterate over
effect_types <- c("shared", "specific", "patterned")
titles <- c("Shared effects", "Trait-specific effects", "Patterned effects")

# Use fixed parameters
noise_type <- "Str. noise"
prop_real_effects <- 0.01
corr_strength <- list(cov1 = 0.3, cov2 = 0.7)

# Modified simulate function to take effect_type as input
simulate_effect_pattern <- function(effect_type) {
  set.seed(123)
  n_genes <- 18000
  n_traits <- 10
  n_real <- round(n_genes * prop_real_effects)
  B_true <- matrix(0, n_genes, n_traits)
  
  # Covariance matrices
  cov_shared <- matrix(0.8, n_traits, n_traits); diag(cov_shared) <- 1
  cov_group1 <- matrix(0.1, n_traits, n_traits)
  cov_group2 <- matrix(0.1, n_traits, n_traits)
  cov_group1[1:5, 1:5] <- corr_strength$cov1
  cov_group2[6:10, 6:10] <- corr_strength$cov2
  diag(cov_group1) <- diag(cov_group2) <- 1
  
  # Populate B_true based on effect type
  if (effect_type == "shared") {
    B_true[1:n_real, ] <- rmvnorm(n_real, rep(0, n_traits), cov_shared)
  } else if (effect_type == "specific") {
    for (i in 1:n_real) {
      trait <- sample(1:n_traits, 1)
      B_true[i, trait] <- rnorm(1)
    }
  } else if (effect_type == "patterned") {
    n_group1 <- round(n_real / 2)
    n_group2 <- n_real - n_group1
    B_true[1:n_group1, ] <- rmvnorm(n_group1, rep(0, n_traits), cov_group1)
    B_true[(n_group1 + 1):n_real, ] <- rmvnorm(n_group2, rep(0, n_traits), cov_group2)
  }
  
  # Generate SE and noise
  se_factors <- runif(n_traits, 0.8, 1.2)
  SE <- matrix(runif(n_genes * n_traits, 0.1, 1), n_genes, n_traits)
  SE <- sweep(SE, 2, se_factors, "*")
  
  # Noise structure: structured unless trait-specific
  if (effect_type == "specific") {
    # Unstructured noise: identity correlation (diagonal)
    R_estimation <- diag(n_traits)
  } else {
    # Structured noise
    random_matrix <- matrix(rnorm(n_traits^2), n_traits)
    Sigma_raw <- crossprod(random_matrix)
    D_inv <- diag(1 / sqrt(diag(Sigma_raw)))
    R_estimation <- D_inv %*% Sigma_raw %*% D_inv
  }
  
  B_obs <- matrix(0, n_genes, n_traits)
  for (i in 1:n_genes) {
    Sigma_i <- diag(SE[i, ]) %*% R_estimation %*% diag(SE[i, ])
    B_obs[i, ] <- rmvnorm(1, B_true[i, ], Sigma_i)
  }
  
  # Run mash
  mash_data <- mash_set_data(B_obs, SE)
  m_1by1 <- mash_1by1(mash_data)
  strong_subset <- get_significant_results(m_1by1, 0.05)
  Vhat <- estimate_null_correlation_simple(mash_data)
  mash_strong <- mash_set_data(B_obs[strong_subset, ], SE[strong_subset, ], V = Vhat)
  U_pca <- cov_pca(mash_strong, 5)
  U_ed <- cov_ed(mash_strong, c(U_pca))
  U_cov <- cov_canonical(mash_data)
  mash_result <- mash(mash_data, Ulist = c(U_ed, U_cov))
  
  # RRMSE
  rrmse <- function(est, truth) sqrt(mean((est - truth)^2)) / sqrt(mean(truth^2))
  B_post <- mash_result$result$PosteriorMean
  rrmse_vals <- c("Single-trait" = rrmse(as.vector(B_obs), as.vector(B_true)),
                  "Multi-trait"  = rrmse(as.vector(B_post), as.vector(B_true)))
  df_rrmse <- data.frame(type = factor(names(rrmse_vals), levels = c("Single-trait", "Multi-trait")),
                         value = rrmse_vals)
  p_rrmse <- ggplot(df_rrmse, aes(x = type, y = value, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(value, 3)), vjust = -0.5, size = 6) +
    labs(x = "Estimate Type", y = "Relative RMSE") +
    scale_fill_manual(values = c("firebrick", "darkblue")) +
    theme_minimal(base_size = 16) +
    scale_y_continuous(limits = c(0, 20)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title.x = element_text(size = 18, color = "black"),
          axis.title.y = element_text(size = 18, color = "black"))
  
  # Lambda plot
  lambda_prior <- get_lambda_singletrait(B_obs, SE)
  lambda_post  <- get_lambda_post(mash_result)
  summary_ci <- function(x) {
    mean_val <- mean(x)
    stderr <- sd(x) / sqrt(length(x))
    c(mean = mean_val, lower = mean_val - 1.96 * stderr, upper = mean_val + 1.96 * stderr)
  }
  df_lambda <- data.frame(
    type = factor(c("Single-trait", "Multi-trait"), levels = c("Single-trait", "Multi-trait")),
    rbind(summary_ci(lambda_prior$Inflation), summary_ci(lambda_post$Inflation))
  )
  p_lambda <- ggplot(df_lambda, aes(x = type, y = mean, fill = type)) +
    geom_hline(aes(yintercept = 1)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_text(aes(label = round(mean, 3)), vjust = -1.5, size = 6) +
    labs(x = "Estimate Type", y = expression("Mean Genomic Inflation Lambda (95% CI)")) +
    scale_fill_manual(values = c("firebrick", "darkblue")) +
    theme_minimal(base_size = 16) +
    scale_y_continuous(limits = c(0, 2)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title.x = element_text(size = 18, color = "black"),
          axis.title.y = element_text(size = 18, color = "black"))
  
  # ROC plot
  truth_vector <- as.vector(abs(B_true) > 1e-6)
  roc_post  <- roc(truth_vector, as.vector(abs(B_post / mash_result$result$PosteriorSD)))
  roc_prior <- roc(truth_vector, as.vector(abs(B_obs / SE)))
  df_post <- data.frame(fpr = 1 - roc_post$specificities, tpr = roc_post$sensitivities)
  df_prior <- data.frame(fpr = 1 - roc_prior$specificities, tpr = roc_prior$sensitivities)
  auc_post <- round(auc(roc_post), 3)
  auc_prior <- round(auc(roc_prior), 3)
  ci_post <- ci.auc(roc_post); ci_prior <- ci.auc(roc_prior)
  legend_labels <- c(
    paste0("Multi-trait \n(AUC = ", auc_post, ", 95% CI:\n", round(ci_post[1], 3), "-", round(ci_post[3], 3), ")"),
    paste0("Single-trait \n(AUC = ", auc_prior, ", 95% CI:\n", round(ci_prior[1], 3), "-", round(ci_prior[3], 3), ")")
  )
  p_roc_main <- ggplot() +
    geom_line(data = df_post, aes(x = fpr, y = tpr, color = "Multi-trait"), size = 1.5) +
    geom_line(data = df_prior, aes(x = fpr, y = tpr, color = "Single-trait"), size = 1.5) +
    geom_abline(linetype = "dashed") +
    labs(x = "False Positive Rate", y = "True Positive Rate", color = NULL) +
    scale_color_manual(values = c("Multi-trait" = "darkblue", "Single-trait" = "firebrick"),
                       labels = legend_labels) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "bottom",
          legend.justification = "left",
          legend.direction = "vertical",
          legend.text = element_text(size = 18),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title.x = element_text(size = 18, color = "black"),
          axis.title.y = element_text(size = 18, color = "black"))
  p_roc <- plot_grid(p_roc_main, nrow = 2, rel_heights = c(1, 0.15))
  
  list(lambda = p_lambda, rrmse = p_rrmse, roc = p_roc)
}

#===========================#
# Run Simulation & Plotting
#===========================#
set.seed(123)
results <- lapply(effect_types, simulate_effect_pattern)

lambda_plots <- lapply(seq_along(results), function(i) {
  ggdraw() + draw_plot(results[[i]]$lambda) +
    draw_label(titles[i], fontface = 'bold', size = 18, x = 0.5, y = 1, hjust = 0.5, vjust = 1)
})
rrmse_plots <- lapply(seq_along(results), function(i) {
  ggdraw() + draw_plot(results[[i]]$rrmse) +
    draw_label(titles[i], fontface = 'bold', size = 18, x = 0.5, y = 1, hjust = 0.5, vjust = 1)
})
roc_plots <- lapply(seq_along(results), function(i) {
  ggdraw() + draw_plot(results[[i]]$roc) +
    draw_label(titles[i], fontface = 'bold', size = 18, x = 0.5, y = 1, hjust = 0.5, vjust = 1)
})

#===========================#
# Output Figures to PDF
#===========================#
pdf("lambda_comparison_effect_types.pdf", width = 18, height = 6)
plot_grid(plotlist = lambda_plots, ncol = 3)
dev.off()

pdf("rrmse_comparison_effect_types.pdf", width = 18, height = 6)
plot_grid(plotlist = rrmse_plots, ncol = 3)
dev.off()

pdf("roc_comparison_effect_types.pdf", width = 18, height = 6)
plot_grid(plotlist = roc_plots, ncol = 3)
dev.off()
