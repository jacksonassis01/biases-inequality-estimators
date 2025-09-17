# Configuration ================================================================
rm(list = ls(all.names = TRUE))

library(mixtools)  
library(ggplot2)   
library(hypergeo)  
library(gtools)    
library(ggplot2)   
library(stats)
library(parallel)
library(tictoc)

source("functions.r")

#===============================================================================

alpha_values <- list(
  c(0.5, 0.5),
  c(0.5, 2.0),
  c(0.5, 3.0),
  c(0.5, 4.0),
  c(0.5, 5.0)
)

n_sim <- 9 # 100  #### Number of Monte Carlo simulations
sample_size <- 15 #20 #### Fixed sample size
true_pi <- c(0.60, 0.40)
true_beta <- c(1, 1)

results <- list()

n_cores = detectCores() - 1

generate_sample_gamma_mixture = function(n, pi, alpha, beta) {
  z <- sample(1:2, size = n, replace = TRUE, prob = pi)
  data_gamma <- rep(NA, n)
  for (i in 1:2) {
    data_gamma[z == i] <- rgamma(sum(z == i), shape = alpha[i], rate = beta[i])
  }
  
  return (data_gamma)
}

tic("parallel")

for (alpha_set in alpha_values) {

  sim_results <- mclapply(1:n_sim, function(sim) {
    data_gamma = generate_sample_gamma_mixture(
      sample_size,
      true_pi,
      alpha_set,
      true_beta
    )
    
    em_out <- estimate_gamma_mixture(data_gamma)
    
    est_pi <- c(em_out$pi1, em_out$pi2)
    est_alpha <- c(em_out$alpha1, em_out$alpha2)
    est_beta <- c(em_out$beta)
    
    est_Gini <- G_estimator(X = data_gamma)
    est_Bias_G <- compute_bias_Ghat_G(n = sample_size, pii = est_pi, alpha = est_alpha)
    est_Corrected_Gini <- est_Gini - est_Bias_G
    
    list(pi = est_pi, alpha = est_alpha, beta = est_beta,
         Gini = est_Gini, Bias_G = est_Bias_G, Corr_Gini = est_Corrected_Gini)
  }, mc.cores = n_cores)
  
  # agora agregamos os resultados
  est_pi <- do.call(rbind, lapply(sim_results, `[[`, "pi"))
  est_alpha <- do.call(rbind, lapply(sim_results, `[[`, "alpha"))
  est_beta <- do.call(rbind, lapply(sim_results, `[[`, "beta"))
  est_Gini <- sapply(sim_results, `[[`, "Gini")
  est_Bias_G <- sapply(sim_results, `[[`, "Bias_G")
  est_Corr_Gini <- sapply(sim_results, `[[`, "Corr_Gini")
  
  results[[paste0("alpha_", paste(alpha_set, collapse = "_"))]] <- list(
    Gini_Mean = mean(est_Gini),
    Bias_G_Mean = mean(est_Bias_G),
    Corrected_Gini_Mean = mean(est_Corr_Gini),
    True_Gini = compute_Gini_G(pii = true_pi, alpha = alpha_set)
  )
}

toc()

# n=1 34
# n=9 180 / 163 / 213 / 169 / 194