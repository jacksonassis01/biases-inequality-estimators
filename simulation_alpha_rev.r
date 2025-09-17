# Configuration ================================================================
rm(list = ls(all.names = TRUE))

library(mixtools)  
library(ggplot2)   
library(hypergeo)  
library(gtools)    
library(ggplot2)   
library(stats)

source("functions.r")

#===============================================================================

alpha_values <- list(
  c(0.5, 0.5),
  c(0.5, 2.0),
  c(0.5, 3.0),
  c(0.5, 4.0),
  c(0.5, 5.0)
)

n_sim <- 100  #### Number of Monte Carlo simulations
sample_size <- 15 #20   #### Fixed sample size
true_pi <- c(0.60, 0.40)
true_beta <- c(1, 1)

results <- list()

################################################################################
##### Monte Carlo Simulation for Different Alpha Values
################################################################################

for (alpha_set in alpha_values) {
  est_pi <- matrix(NA, nrow = n_sim, ncol = 2)
  est_alpha <- matrix(NA, nrow = n_sim, ncol = 2)
  est_beta <- matrix(NA, nrow = n_sim, ncol = 2)
  est_Gini <- rep(NA, n_sim)
  est_Bias_G <- rep(NA, n_sim)
  est_Corrected_Gini <- rep(NA, n_sim)
  
  for (sim in 1:n_sim) {
    data_gamma = generate_sample_gamma_mixture(
      sample_size,
      true_pi,
      alpha_set,
      true_beta
    )
    
    em_out <- estimate_gamma_mixture(data_gamma)
    
    est_pi[sim, ] <-  c(em_out$pi1, em_out$pi2)
    est_alpha[sim, ] <- c(em_out$alpha1, em_out$alpha2)
    est_beta[sim, ] <-  c(em_out$beta)
    
    est_Gini[sim] <- G_estimator(X = data_gamma)
    
    est_Bias_G[sim] <- compute_bias_Ghat_G(
      n = sample_size,
      pii = est_pi[sim, ],
      alpha = est_alpha[sim, ]
    )
    
    est_Corrected_Gini[sim] <- est_Gini[sim] - est_Bias_G[sim]
  }
  
  results[[paste0("alpha_", paste(alpha_set, collapse = "_"))]] <- list(
    Gini_Mean = mean(est_Gini),
    Bias_G_Mean = mean(est_Bias_G),
    Corrected_Gini_Mean = mean(est_Corrected_Gini),
    True_Gini = compute_Gini_G(pii = true_pi, alpha = alpha_set)
  )
}

################################################################################
##### Print Results
################################################################################

for (alpha_label in names(results)) {
  cat("\nResults for Alpha Set:", alpha_label, "\n")
  cat("Mean Gini Estimate:", results[[alpha_label]]$Gini_Mean, "\n")
  cat("Mean Bias in Gini Estimate:", results[[alpha_label]]$Bias_G_Mean, "\n")
  cat("Mean Corrected Gini Estimate:", results[[alpha_label]]$Corrected_Gini_Mean, "\n")
  cat("True Gini:", results[[alpha_label]]$True_Gini, "\n")
}

# Generate Plot ================================================================
data_results <- data.frame(
  Alpha_Set = c(1, 2, 3, 5), #c(0.5, 1, 2, 3, 5),
  Mean_Gini_Estimate = sapply(results, function(x) x$Gini_Mean),
  Mean_Corrected_Gini = sapply(results, function(x) x$Corrected_Gini_Mean)
)

data_long <- reshape2::melt(data_results, id.vars = "Alpha_Set")

#### Convertendo para formato longo para facilitar o gráfico

#### Gerando o gráfico
plot <- ggplot(data_long, aes(x = Alpha_Set, y = value, color = variable)) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = expression(alpha), y = "average value", 
       title = "") +
  theme_minimal() +
  scale_color_manual(name = "Estimator", values = c("red","blue"),
                     labels = c("Gini", "Bias-corrected Gini")) +
   theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

#### Exibir o gráfico
print(plot)

ggsave(
  filename = "gini_comparison_alphas.eps",
  plot = plot,
  device = cairo_ps,
  width = 7,
  height = 5
)

#===============================================================================
