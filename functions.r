# Function to compute the Gini coefficient G ===================================
G_estimator <- function(X) {
  n <- length(X)
  if (n < 2) {
    stop("Sample size must be at least 2.")
  }
  
  numerator <- sum(outer(X, X, function(xi, xj) abs(xi - xj))) / 2
  denominator <- sum(X)
  
  G_hat <- (1 / (n - 1)) * (numerator / denominator)
  
  return(G_hat)
}

compute_Gini_G <- function(pii, alpha) {
  
  m <- length(pii)  #### Number of components
  
  #### Compute π*
  pi_star <- pii * alpha / sum(pii * alpha)
  
  #### Compute Beta function values
  B_ij <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      B_ij[i, j] <- beta(alpha[i], alpha[j] + 1)
    }
  }
  
  #### Compute 2F1 hypergeometric function at 1/2
  hypergeom_values <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      hypergeom_values[i, j] <- as.numeric( hypergeo(alpha[i], -alpha[j], alpha[i] + 1, 1/2))
    }
  }
  
  #### Compute Gini coefficient G
  G <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      G <- G + (pii[i] * pi_star[j] * hypergeom_values[i, j] / (2^(alpha[i] - 1) * alpha[i] * B_ij[i, j]) )
    }
  }
  
  G <- G - 1  #### Apply the final adjustment in Equation (7)
  
  return(G)
}

#===============================================================================

# Function to compute Bias(Ĝ, G) ===============================================
compute_bias_Ghat_G <- function(n, pii, alpha) {
  
  m <- length(pii)  #### Number of components
  
  #### Compute π*
  pi_star <- pii * alpha / sum(pii * alpha)
  
  #### Compute Beta function values
  B_ij <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      B_ij[i, j] <- beta(alpha[i], alpha[j] + 1)
    }
  }
  
  #### Compute 2F1 hypergeometric function at 1/2
  hypergeom_values <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      hypergeom_values[i, j] <- as.numeric(hypergeo(-alpha[j], alpha[i], alpha[i] + 1, 1/2))
    }
  }
  
  #### Compute First Bias Term
  Bias1_G_hat_G <- 0
  
  for (i in 1:m) {
    for (j in 1:m) {
      term1 <- 0
      
      #### Generate all index combinations for (n-2) elements
      subset_indices <- expand.grid(rep(list(1:m), (n-2)))
      
      for (k in 1:nrow(subset_indices)) {
        subset_vector <- as.numeric(subset_indices[k, ])  #### Convert to numeric vector
        subset_alpha  <- sum(alpha[subset_vector])
        subset_pi_prod <- prod(pii[subset_vector])
        term1 <- term1 + ((subset_pi_prod * n * sum(pii * alpha)) / (subset_alpha + alpha[i] + alpha[j]))    
      }
      
      term1 <- term1 - 1 
      
      Bias1_G_hat_G <- Bias1_G_hat_G + (pii[i] * pi_star[j] * hypergeom_values[i, j] / 
                                          (2^(alpha[i] - 1) * alpha[i] * B_ij[i, j])) * term1
    }
  }
  
  #### Compute Second Bias Term
  Bias2_G_hat_G <- 0
  
  #### Generate all index combinations for (n-1) elements
  subset_indices <- expand.grid(rep(list(1:m), (n-1)))
  
  for (k in 1:nrow(subset_indices)) {
    subset_vector <- as.numeric(subset_indices[k, ])  #### Convert to numeric vector
    subset_alpha <- sum(alpha[subset_vector])
    subset_pi_prod <- prod(pii[subset_vector])
    term2 <- 0
    for (j in 1:m) {
      term2 <- term2 + (n * pii[j] * alpha[j]) / (subset_alpha + alpha[j])
    }
    Bias2_G_hat_G <- Bias2_G_hat_G + subset_pi_prod * term2
  }
  
  Bias2_G_hat_G <- 1 - Bias2_G_hat_G
  
  #### Final Bias Value
  Bias_G_hat_G <- Bias1_G_hat_G + Bias2_G_hat_G
  
  return(Bias_G_hat_G)
}

bias_corrected_Ghat <- function(X, n, pii, alpha) {
  m <- length(pii)  #### Number of components
  n <- length(X)
  
  corrected_Ghat_value <- G_estimator(X = X) - compute_bias_Ghat_G(n = n, pii = pii , alpha =alpha)
  
  return(corrected_Ghat_value)
}

#===============================================================================

# Function to estimate gamma mixture parameters ================================

# Função de verossimilhança negativa
neg_log_likelihood <- function(params, data) {
  #### Parâmetros a serem estimados
  pi1 <- params[1]  #### Proporção da primeira componente
  alpha1 <- params[2]  #### Parâmetro shape da primeira Gamma
  alpha2 <- params[3]  #### Parâmetro shape da segunda Gamma
  beta <- params[4]  #### Parâmetro rate (comum para ambas as componentes)
  
  #### Garante que pi1 esteja no intervalo [0,1]
  pi1 <- exp(pi1) / (1 + exp(pi1))  
  pi2 <- 1 - pi1  #### Como são duas misturas, pi2 = 1 - pi1
  
  #### Função de densidade da mistura de duas Gammas
  density_mixture <- pi1 * dgamma(data, shape = alpha1, rate = beta) + 
    pi2 * dgamma(data, shape = alpha2, rate = beta)
  
  #### Evita log(0) usando um limite inferior muito pequeno
  density_mixture[density_mixture < 1e-10] <- 1e-10
  
  #### Retornar o valor negativo da soma dos log-verossimilhanças
  return(-sum(log(density_mixture)))
}

#### Função para estimar os parâmetros via MLE
estimate_gamma_mixture <- function(data, start_params = c(0.5, 2, 5, 1)) {
  #### Otimização via optim()
  mle_fit <- optim(
    par = c(log(start_params[1] / (1 - start_params[1])), start_params[2], start_params[3], start_params[4]),  #### Transformação logit para pi1
    fn = neg_log_likelihood,
    data = data,
    method = "L-BFGS-B",
    lower = c(-Inf, 0.001, 0.001, 0.001),  #### Limita alpha > 0 e beta > 0
    upper = c(Inf, Inf, Inf, Inf)
  )
  
  #### Recuperando os parâmetros estimados
  pi1_est <- exp(mle_fit$par[1]) / (1 + exp(mle_fit$par[1]))  #### Inverso da transformação logit
  pi2_est <- 1 - pi1_est
  alpha1_est <- mle_fit$par[2]
  alpha2_est <- mle_fit$par[3]
  beta_est <- mle_fit$par[4]
  
  #### Retornar os parâmetros estimados
  return(list(
    pi1 = pi1_est, 
    pi2 = pi2_est, 
    alpha1 = alpha1_est, 
    alpha2 = alpha2_est, 
    beta = beta_est, 
    logLik = -mle_fit$value
  ))
}

#===============================================================================

generate_sample_gamma_mixture = function(n, pi, alpha, beta) {
  z <- sample(1:2, size = n, replace = TRUE, prob = pi)
  data_gamma <- rep(NA, n)
  for (i in 1:2) {
    data_gamma[z == i] <- rgamma(sum(z == i), shape = alpha[i], rate = beta[i])
  }
  
  return (data_gamma)
}