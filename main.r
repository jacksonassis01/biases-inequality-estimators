# configuration ================================================================
rm(list = ls())

library(gtools)

source("functions.r")

SEED = 123

#===============================================================================

# parameters ===================================================================
pis <- c(.5, .5)
alphas <- c(1, 2)
beta <- 1
ns <- c(10, 100, 1000, 10000)

#===============================================================================

# [x] funcao para gerar os dados ===============================================
generate_sample_gamma_mixture <- function(n, pi, alpha, beta, seed = 123) {
  set.seed(seed)
  z <- sample(1:2, size = n, replace = TRUE, prob = pi)
  data_gamma <- rep(NA, n)
  for (i in 1:2) {
    data_gamma[z == i] <- rgamma(sum(z == i), shape = alpha[i], rate = beta[i])
  }
  
  return (data_gamma)
}

z <- generate_sample_gamma_mixture(30, c(.5, .5), c(1, 4), c(1, 1))

hist(z)

par(mfrow = c(1,3))

set.seed(SEED)
x <- rgamma(30, 1, 1)
y <- rgamma(30, 1, 4)

hist(z)
hist(x)
hist(y)

par(mfrow = c(1,1))

#===============================================================================

# [x] funcao para estimar os parametros ========================================
estimate_parameters <- function(data, start_params = c(0.5, 2, 5, 1)) {
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

fit <- estimate_parameters(z)

fits <- lapply(ns, \(n) {
  z <- generate_sample_gamma_mixture(n, c(.5, .5), c(1, 4), c(1, 1))
  fit <- estimate_parameters(z)
  return(fit)
})

cbind(ns, do.call(rbind, fits))

#===============================================================================

# [x] funcoes dos estimadores ==================================================
Thiel_T_estimator <- function(X) {
  mu <- mean(X)
  T_T <- mean((X / mu) * log(X / mu))
  return(T_T)
}

TL_estimator <- function(X) {
  mu <- mean(X)
  T_L <- mean(log(mu / X))
  return(T_L)
}

A1_estimator <- function(X) {
  mu <- mean(X)
  gm <- exp(mean(log(X)))
  return(1 - gm / mu)
}

AInf_estimator <- function(X) {
  mu <- mean(X)
  return(1 - min(X) / mu)
}

DI_estimator <- function(X) {
  mu <- mean(X)
  return(var(X) / mu)
}

estimates <- lapply(ns, \(n) {
  z <- generate_sample_gamma_mixture(n, c(.5, .5), c(1, 4), c(1, 1))
  
  return(list(
    Thiel_T = Thiel_T_estimator(z),
    TL = TL_estimator(z),
    A1 = A1_estimator(z),
    AInf = AInf_estimator(z),
    DI = DI_estimator(z)
  ))
})

cbind(ns, do.call(rbind, estimates))

#===============================================================================

# [x] computar os valores verdadeiros dos indices para vies / eqm ==============
compute_Theil_T <- function(pi, alpha, lambda = 1) {
  S <- sum(pi * alpha)
  sum_pi_alpha_psi <- sum(pi * alpha * digamma(alpha))
  Thiel_T <- (1 / S) * (sum_pi_alpha_psi + 1) - log(S)
  return(Thiel_T)
}

compute_Theil_L <- function(pi, alpha, lambda = 1) {
  S <- sum(pi * alpha)
  sum_pi_psi <- sum(pi * digamma(alpha))
  TL <- log(S) - sum_pi_psi
  return(TL)
}

compute_Atkinson_e1 <- function(pi, alpha, lambda = 1) {
  S <- sum(pi * alpha)
  sum_pi_psi <- sum(pi * digamma(alpha))
  A1 <- 1 - (1 / S) * exp(sum_pi_psi)
  return(A1)
}

compute_Atkinson_eInf <- function(pi, alpha, lambda = 1) {
  return(1)
}

compute_DI <- function(pi, alpha, lambda = 1) {
  S <- sum(pi * alpha)
  sum_pi_alpha_alpha1 <- sum(pi * alpha * (alpha + 1))
  VMR <- (1 / (lambda * S)) * (sum_pi_alpha_alpha1 - S^2)
  return(VMR)
}

compute_true_values_indices = function(pis, alphas, beta) {
  return(list(
    Thiel_T = compute_Theil_T(pis, alphas, beta),
    TL = compute_Theil_L(pis, alphas, beta),
    A1 = compute_Atkinson_e1(pis, alphas, beta),
    AInf = compute_Atkinson_eInf(pis, alphas, beta),
    DI = compute_DI(pis, alphas, beta)
  ))
}

true_values_indices = compute_true_values_indices(pis, alphas, beta)

errors = lapply(estimates, \(estimates_for_n) {
  sapply(names(estimates_for_n), function(name) {
    estimates_for_n[[name]] - true_values_indices[[name]]
  })
})

do.call(rbind, errors)

#===============================================================================

# [x] computar os valores dos vieses ===========================================
compute_bias_Thiel_That <- function(pi, alpha, n) {
  m <- length(pi)
  
  indices <- permutations(n = m, r = n, v = 1:m, repeats.allowed = TRUE)
  
  sum1 <- 0
  for (k in 1:nrow(indices)) {
    j <- indices[k, ]
    pi_prod <- prod(pi[j])
    alpha_sum <- sum(alpha[j])
    sum1 <- sum1 + pi_prod * ( (sum(alpha[j] * digamma(alpha[j])) - lgamma(alpha_sum) + n - 1) / alpha_sum )
  }
  
  sum1 <- sum1 + log(n)
  
  alpha_bar <- sum(pi * alpha)
  sum2 <- (sum(pi * alpha * digamma(alpha)) + 1) / alpha_bar - log(alpha_bar)
  
  bias <- sum1 - sum2
  
  return(bias)
}

#===============================================================================

# [ ] funcao para comparar o desempenho do classico com o do vies ==============
compare_estimators_mixture <- function(pi, alpha, beta, n, estimators, true_values, bias_funcs = NULL, n_sim = 1000) {
  results <- list()
  
  for (name in names(estimators)) {
    est_classic <- numeric(n_sim)
    est_biascorr <- numeric(n_sim)
    
    for (s in 1:n_sim) {
      # 1. Gerar amostra da mistura de gamas
      X <- generate_sample_gamma_mixture(n, pi, alpha, beta)
      
      # 2. Estimador clássico
      est_classic[s] <- estimators[[name]](X)
      
      # 3. Estimador com correção de viés (se fornecida)
      if (!is.null(bias_funcs) && !is.null(bias_funcs[[name]])) {
        bias_val <- bias_funcs[[name]](pi, alpha, n)
        est_biascorr[s] <- est_classic[s] - bias_val
      } else {
        est_biascorr[s] <- NA
      }
    }
    
    # Valor verdadeiro fornecido
    true_val <- true_values[[name]]
    
    # Armazenar resultados
    results[[name]] <- list(
      classic = list(
        bias = mean(est_classic, na.rm = TRUE) - true_val,
        var  = var(est_classic, na.rm = TRUE),
        mse  = mean((est_classic - true_val)^2, na.rm = TRUE)
      ),
      bias_corrected = list(
        bias = mean(est_biascorr, na.rm = TRUE) - true_val,
        var  = var(est_biascorr, na.rm = TRUE),
        mse  = mean((est_biascorr - true_val)^2, na.rm = TRUE)
      )
    )
  }
  
  return(results)
}

estimators = list(
  Thiel_T = Thiel_T_estimator
)

bias_funcs = list(
  Thiel_T = compute_bias_Thiel_That
  
)

compare_estimators_mixture(
  pi = pis,
  alpha = alphas,
  beta = c(1, 1),
  n = 15,
  estimators = estimators,
  true_values = true_values_indices,
  bias_funcs = bias_funcs,
  n_sim = 10
)

#===============================================================================

# [ ] contrastar os estimadores naturais e os corridos =========================
#===============================================================================

# [ ] variar os cenarios (m, peso da mistura igual, alphas) ====================
#===============================================================================

# [ ] n = 20 ===================================================================
#===============================================================================

# [ ] comparar a performance com e sem correcao de vies ========================
#===============================================================================
