########################################################################
# Title: Simulations to Evaluate mash Performance Under Varying Settings
# Description: 
#   This script simulates multivariate genetic association data under 
#   various noise structures, effect sparsity levels, and trait 
#   correlation patterns. It evaluates performance of single-trait vs. 
#   mash (multi-trait) analysis using genomic inflation (lambda), 
#   relative RMSE, and ROC-AUC.
# Author: Kiran Biddinger
# Date: April 2025
########################################################################

#=======================#
# Load Required Libraries
#=======================#
library(ashr)
library(mashr)
library(mvtnorm)
library(ggplot2)
library(gaston)
library(MASS)
library(pROC)
library(gridExtra)
library(grid)
library(cowplot)

#=======================#
# Set Output Directory
#=======================#
setwd("/Users/kbidding/Documents/mashHF/Figures/Simulations")

#=======================#
# Simulation Parameters
#=======================#
noise_types   <- c("Str. noise", "Unstr. noise")
prop_effects  <- c(0.005, 0.01)
corr_types    <- list(
  "strong trait corr." = list(cov1 = 0.9, cov2 = 0.7),
  "weak trait corr."   = list(cov1 = 0.3, cov2 = 0.1)
)

#====================================================#
# Helper Function: Compute Lambda from Observed Data
#====================================================#
get_lambda_singletrait <- function(B_obs, SE) {
  n_col <- ncol(B_obs)
  results <- data.frame()
  
  for (i in 1:n_col) {
    beta = as.data.frame(B_obs)[, i]
    se = as.data.frame(SE)[, i]
    z    <- beta / se
    chisq <- z^2
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    lambda_95 <- quantile(chisq, 0.95, na.rm = TRUE) / qchisq(0.05, df = 1, lower.tail = FALSE)
    
    results <- rbind(results, data.frame(Trait = paste0("Trait_", i), 
                                         Inflation = lambda_50, 
                                         Inflation_95 = lambda_95))
  }
  
  colnames(results) <- c("Trait", "Inflation", "Inflation_95")
  results$Inflation <- as.numeric(results$Inflation)
  return(results)
}

#====================================================#
# Helper Function: Compute Lambda from mash Posterior
#====================================================#
get_lambda_post <- function(data) {
  n_col <- ncol(data$result$PosteriorMean)
  colnames(data$result$PosteriorMean) <- paste0("BETA_", seq_len(n_col))
  colnames(data$result$PosteriorSD)   <- paste0("SE_", seq_len(n_col))
  data$result$PosteriorMean <- as.data.frame(data$result$PosteriorMean)
  data$result$PosteriorSD <- as.data.frame(data$result$PosteriorSD)
  
  results <- data.frame()
  for (i in 1:n_col) {
    beta <- data$result$PosteriorMean[[paste0("BETA_", i)]]
    se   <- data$result$PosteriorSD[[paste0("SE_", i)]]
    z    <- beta / se
    chisq <- z^2
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    results <- rbind(results, data.frame(Trait = paste0("Trait_", i), Inflation = lambda_50))
  }
  
  colnames(results) <- c("Trait", "Inflation")
  results$Inflation <- as.numeric(results$Inflation)
  return(results)
}

#====================================================#
# Function: Simulate a Setting and Generate 3 Plots
#====================================================#
simulate_and_plot <- function(noise_type, prop_real_effects, corr_strength) {
  n_genes <- 18000
  n_traits <- 10
  
  # Generate covariance matrices
  cov_shared <- matrix(0.8, n_traits, n_traits); diag(cov_shared) <- 1
  cov_group1 <- matrix(0.1, n_traits, n_traits)
  cov_group2 <- matrix(0.1, n_traits, n_traits)
  cov_group1[1:5, 1:5] <- corr_strength$cov1
  cov_group2[6:10, 6:10] <- corr_strength$cov2
  diag(cov_group1) <- diag(cov_group2) <- 1
  
  # Simulate true effects
  n_real <- round(n_genes * prop_real_effects)
  n_shared <- round(n_real * 0.3)
  n_specific <- round(n_real * 0.4)
  n_group1 <- round(n_real * 0.15)
  n_group2 <- n_real - (n_shared + n_specific + n_group1)
  
  B_true <- matrix(0, n_genes, n_traits)
  if (n_shared > 0) B_true[1:n_shared, ] <- rmvnorm(n_shared, rep(0, n_traits), cov_shared)
  if (n_specific > 0) {
    for (i in (n_shared + 1):(n_shared + n_specific)) {
      trait <- sample(1:n_traits, 1)
      B_true[i, trait] <- rnorm(1)
    }
  }
  if (n_group1 > 0) B_true[(n_shared + n_specific + 1):(n_shared + n_specific + n_group1), ] <- rmvnorm(n_group1, rep(0, n_traits), cov_group1)
  if (n_group2 > 0) B_true[(n_shared + n_specific + n_group1 + 1):(n_real), ] <- rmvnorm(n_group2, rep(0, n_traits), cov_group2)
  
  # Simulate standard errors
  se_factors <- runif(n_traits, 0.8, 1.2)
  SE <- matrix(runif(n_genes * n_traits, 0.1, 1), n_genes, n_traits)
  SE <- sweep(SE, 2, se_factors, "*")
  
  # Generate noise structure
  R_estimation <- if (tolower(noise_type) == "str. noise") {
    Sigma_raw <- crossprod(matrix(rnorm(n_traits^2), n_traits))
    D_inv <- diag(1 / sqrt(diag(Sigma_raw)))
    D_inv %*% Sigma_raw %*% D_inv
  } else diag(n_traits)
  
  # Simulate observed effects
  B_obs <- matrix(0, n_genes, n_traits)
  for (i in 1:n_genes) {
    Sigma_i <- diag(SE[i, ]) %*% R_estimation %*% diag(SE[i, ])
    B_obs[i, ] <- rmvnorm(1, B_true[i, ], Sigma_i)
  }
  
  # Fit mash
  mash_data <- mash_set_data(B_obs, SE)
  m_1by1 <- mash_1by1(mash_data)
  strong_subset <- get_significant_results(m_1by1, 0.1)
  Vhat <- estimate_null_correlation_simple(mash_data)
  mash_strong <- mash_set_data(B_obs[strong_subset, ], SE[strong_subset, ], V = Vhat)
  U_pca <- cov_pca(mash_strong, 5)
  U_ed <- cov_ed(mash_strong, c(U_pca))
  U_cov <- cov_canonical(mash_data)
  mash_result <- mash(mash_data, Ulist = c(U_ed, U_cov))
  
  # Plot: RRMSE
  rrmse <- function(est, truth) sqrt(mean((est - truth)^2)) / sqrt(mean(truth^2))
  rrmse_vals <- c("Single-trait" = rrmse(as.vector(B_obs), as.vector(B_true)),
                  "Multi-trait"  = rrmse(as.vector(mash_result$result$PosteriorMean), as.vector(B_true)))
  df_rrmse <- data.frame(type = factor(names(rrmse_vals), levels = c("Single-trait", "Multi-trait")), 
                         value = rrmse_vals)
  p_rrmse <- ggplot(df_rrmse, aes(x = type, y = value, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(value, 3)), vjust = -0.5, size = 6) +
    labs(x = NULL, y = "Relative RMSE") +
    scale_fill_manual(values = c("firebrick", "darkblue")) +
    theme_minimal(base_size = 18) +
    scale_y_continuous(limits = c(0, 20)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title.x = element_text(size = 18, color = "black"),
          axis.title.y = element_text(size = 18, color = "black"))
  print(p_rrmse)
  
  # Plot: Genomic Inflation (lambda)
  lambda_prior <- get_lambda_singletrait(B_obs, SE)
  lambda_post  <- get_lambda_post(mash_result)
  summary_ci <- function(x) {
    mean_val <- mean(x)
    stderr   <- sd(x) / sqrt(length(x))
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
    labs(x = NULL, y = expression("Mean Genomic Inflation Lambda (95% CI)")) +
    scale_fill_manual(values = c("firebrick", "darkblue")) +
    theme_minimal(base_size = 16) +
    scale_y_continuous(limits = c(0, 2)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title.x = element_text(size = 18, color = "black"),
          axis.title.y = element_text(size = 18, color = "black"))
  print(p_lambda)
  
  # Plot: ROC curves
  truth_vector <- as.vector(abs(B_true) > 1e-6)
  roc_post  <- roc(truth_vector, as.vector(abs(mash_result$result$PosteriorMean / mash_result$result$PosteriorSD)))
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
    scale_color_manual(values = c("Multi-trait" = "darkblue", "Single-trait" = "firebrick"), labels = legend_labels) +
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
  print(p_roc)
  
  list(lambda = p_lambda, roc = p_roc, rrmse = p_rrmse)
}

#=============================#
# Run All Simulation Settings
#=============================#
set.seed(123)
all_plots <- list(lambda = list(), roc = list(), rrmse = list())
setting_labels <- c()
i <- 1
for (noise in noise_types) {
  for (prop in prop_effects) {
    for (corr_name in names(corr_types)) {
      label <- paste0(noise, " and ", corr_name, " (", prop * 100, "% genes true effects)")
      setting_labels[i] <- label
      print(label)
      plots <- simulate_and_plot(noise, prop, corr_types[[corr_name]])
      all_plots$lambda[[i]] <- plots$lambda
      all_plots$roc[[i]] <- plots$roc
      all_plots$rrmse[[i]] <- plots$rrmse
      i <- i + 1
    }
  }
}

#=============================#
# Annotate and Save All Plots
#=============================#
add_title <- function(plot, title) {
  ggdraw() + draw_plot(plot, y = 0.05, height = 0.8) + draw_label(title, fontface = 'bold', size = 18,
                                                                  x = 0.5, y = 0.99, hjust = 0.5, vjust = 1)
}
format_title <- function(label) gsub("\\(", "\n(", label)

plot_titles <- lapply(setting_labels, format_title)
lambda_plots <- mapply(add_title, all_plots$lambda, plot_titles, SIMPLIFY = FALSE)
rrmse_plots  <- mapply(add_title, all_plots$rrmse,  plot_titles, SIMPLIFY = FALSE)
roc_plots    <- mapply(add_title, all_plots$roc,    plot_titles, SIMPLIFY = FALSE)

pdf("lambda_panel_4x2.pdf", width = 16, height = 14)
plot_grid(plotlist = lambda_plots, ncol = 4)
dev.off()

pdf("rrmse_panel_4x2.pdf", width = 16, height = 14)
plot_grid(plotlist = rrmse_plots, ncol = 4, scale = 0.95)
dev.off()

pdf("roc_panel_4x2.pdf", width = 16, height = 14)
plot_grid(plotlist = roc_plots, ncol = 4)
dev.off()
