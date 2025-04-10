########################################################################
# Simulation: Compare Effect Structures Under Biased Noise
# Conditions:
# - (1) Inverse-rank normalized Z-scores
# - (2) Raw Z-scores
# - (3) Full input: beta and standard error
########################################################################

#=======================#
# Load Required Libraries
#=======================#
library(ashr)
library(mashr)
library(mvtnorm)
library(ggplot2)
library(pROC)
library(gridExtra)
library(grid)
library(cowplot)

#=======================#
# Simulation Parameters
#=======================#
set.seed(123)
setwd("/Users/kbidding/Documents/mashHF/Figures/Simulations")
settings <- c("Inverse-rank normalized z-scores", "Raw z-scores", "Full input: Beta and std. error")

#====================================================#
# Helper: Inverse-rank normalization
#====================================================#
irn_transform <- function(z) {
  apply(z, 2, function(trait) {
    qnorm((rank(trait, ties.method = "average") - 0.5) / length(trait))
  })
}

#====================================================#
# Helper: Compute Lambda from Z-score matrix
#====================================================#
get_lambda_z <- function(z) {
  n_col <- ncol(z)
  results <- data.frame()
  
  for (i in 1:n_col) {
    trait_z <- z[, i]
    chisq <- trait_z^2
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    lambda_95 <- quantile(chisq, 0.95, na.rm = TRUE) / qchisq(0.05, df = 1, lower.tail = FALSE)
    results <- rbind(results, data.frame(Trait = paste0("Trait_", i), 
                                         Inflation = lambda_50, 
                                         Inflation_95 = lambda_95))
  }
  return(results)
}

#====================================================#
# Helper: Compute Lambda from beta/se
#====================================================#
get_lambda_singletrait <- function(B_obs, SE) {
  n_col <- ncol(B_obs)
  results <- data.frame()
  
  for (i in 1:n_col) {
    z <- B_obs[, i] / SE[, i]
    chisq <- z^2
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    lambda_95 <- quantile(chisq, 0.95, na.rm = TRUE) / qchisq(0.05, df = 1, lower.tail = FALSE)
    results <- rbind(results, data.frame(Trait = paste0("Trait_", i), 
                                         Inflation = lambda_50, 
                                         Inflation_95 = lambda_95))
  }
  return(results)
}

#====================================================#
# Helper: Compute Lambda from mash posterior
#====================================================#
get_lambda_post <- function(data) {
  post_mean <- data$result$PosteriorMean
  post_sd <- data$result$PosteriorSD
  n_col <- ncol(post_mean)
  
  results <- data.frame()
  for (i in 1:n_col) {
    z <- post_mean[, i] / post_sd[, i]
    chisq <- z^2
    lambda_50 <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
    results <- rbind(results, data.frame(Trait = paste0("Trait_", i), Inflation = lambda_50))
  }
  return(results)
}

#====================================================#
# Main Simulation Function
#====================================================#
simulate_and_plot <- function(setting) {
  n_genes <- 18000
  n_traits <- 10
  prop_real_effects <- 0.01
  corr_strength <- list(cov1 = 0.3, cov2 = 0.7)
  
  #--------------------------#
  # 1. Simulate true effects
  #--------------------------#
  cov_shared <- matrix(0.8, n_traits, n_traits); diag(cov_shared) <- 1
  cov_group1 <- matrix(0.1, n_traits, n_traits)
  cov_group2 <- matrix(0.1, n_traits, n_traits)
  cov_group1[1:5, 1:5] <- corr_strength$cov1
  cov_group2[6:10, 6:10] <- corr_strength$cov2
  diag(cov_group1) <- diag(cov_group2) <- 1
  
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
  
  #--------------------------#
  # 2. Simulate biased SEs and Z-score artifacts
  #--------------------------#
  SE <- matrix(rlnorm(n_genes * n_traits, meanlog = -0.2, sdlog = 0.5), n_genes, n_traits)
  z_bias <- matrix(rnorm(n_genes * n_traits, mean = 0.2, sd = 0.05), n_genes, n_traits)
  B_obs <- B_true + z_bias * SE
  
  #--------------------------#
  # 3. Prepare mash data
  #--------------------------#
  if (setting == "Inverse-rank normalized z-scores") {
    Z <- irn_transform(B_obs / SE)
    mash_data <- mash_set_data(Z)
  } else if (setting == "Raw z-scores") {
    Z <- B_obs / SE
    mash_data <- mash_set_data(Z)
  } else {
    mash_data <- mash_set_data(B_obs, SE)
  }
  if (setting == "Raw z-scores" || setting == "Inverse-rank normalized z-scores") {
    hist(Z, breaks = 100, main = setting, xlab = "Z", col = "skyblue")
  }
  
  #--------------------------#
  # 4. Fit mash
  #--------------------------#
  m_1by1 <- mash_1by1(mash_data)
  strong_subset <- get_significant_results(m_1by1, 0.1)
  if (is.null(strong_subset) || length(strong_subset) == 0) {
    warning("No significant results found. Sampling 100 random genes as fallback.")
    strong_subset <- sample(1:n_genes, 100)
  }  
  if (setting == "Inverse-rank normalized z-scores") {
    mash_strong <- mash_set_data(Z[strong_subset, ])
  } else if (setting == "Raw z-scores") {
    Z <- B_obs / SE
    mash_strong <- mash_set_data(Z[strong_subset, ])
  } else {
    Vhat <- estimate_null_correlation_simple(mash_data)
    mash_strong <- mash_set_data(B_obs[strong_subset, ], SE[strong_subset, ], V=Vhat)
  }  
  U_pca <- cov_pca(mash_strong, 5)
  U_ed <- cov_ed(mash_strong, U_pca)
  U_cov <- cov_canonical(mash_data)
  mash_result <- mash(mash_data, Ulist = c(U_ed, U_cov))
  
  #--------------------------#
  # 5. Lambda plot
  #--------------------------#
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
    geom_hline(yintercept = 1) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_text(aes(label = round(mean, 3)), vjust = -1.5, size = 6) +
    labs(x = "Estimate Type", y = expression("Mean Genomic Inflation Lambda (95% CI)")) +
    scale_fill_manual(values = c("firebrick", "darkblue")) +
    theme_minimal(base_size = 16) +
    ylim(0, 2) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title.x = element_text(size = 18, color = "black"),
          axis.title.y = element_text(size = 18, color = "black"))
  
  #--------------------------#
  # 6. ROC Curve
  #--------------------------#
  truth_vector <- as.vector(abs(B_true) > 1e-6)
  roc_post <- roc(truth_vector, as.vector(abs(mash_result$result$PosteriorMean / mash_result$result$PosteriorSD)))
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
  
  return(list(lambda = p_lambda, roc = p_roc))
}



#=============================#
# Run All Simulation Settings
#=============================#
set.seed(123)
all_plots <- list(lambda = list(), roc = list(), rrmse = list())
for (i in seq_along(settings)) {
  setting <- settings[i]
  print(setting)
  plots <- simulate_and_plot(setting)
  all_plots$lambda[[i]] <- plots$lambda
  all_plots$roc[[i]] <- plots$roc
}
setting_labels <- settings

#=============================#
# Annotate and Save All Plots
#=============================#
add_title <- function(plot, title) {
  ggdraw() + draw_plot(plot, y = 0.05, height = 0.8) + draw_label(title, fontface = 'bold', size = 18,
                                                                  x = 0.5, y = 0.99, hjust = 0.5, vjust = 1)
}

plot_titles <- setting_labels
lambda_plots <- mapply(add_title, all_plots$lambda, plot_titles, SIMPLIFY = FALSE)
rrmse_plots  <- mapply(add_title, all_plots$rrmse,  plot_titles, SIMPLIFY = FALSE)
roc_plots    <- mapply(add_title, all_plots$roc,    plot_titles, SIMPLIFY = FALSE)

pdf("bias_lambda_panel.pdf", width = 16, height = 14)
plot_grid(plotlist = lambda_plots, ncol = 3)
dev.off()

pdf("bias_roc_panel.pdf", width = 16, height = 14)
plot_grid(plotlist = roc_plots, ncol = 3)
dev.off()

