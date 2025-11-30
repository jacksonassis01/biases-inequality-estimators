# R/mixture/mixture_loglik.R
# Verossimilhança para mistura gamma (duas parametrizações)

# Neg log-likelihood (parametrização direta)
neg_log_likelihood <- function(params, data) {
  # params: (pi1_raw, alpha1, alpha2, beta)
  pi1_raw <- params[1]
  alpha1  <- params[2]
  alpha2  <- params[3]
  beta    <- params[4]

  pi1 <- exp(pi1_raw) / (1 + exp(pi1_raw))
  pi2 <- 1 - pi1

  f1 <- dgamma(data, shape = alpha1, rate = beta)
  f2 <- dgamma(data, shape = alpha2, rate = beta)
  density_mixture <- pi1 * f1 + pi2 * f2

  # estabiliza
  density_mixture[density_mixture < 1e-12] <- 1e-12
  -sum(log(density_mixture))
}

# Neg log-likelihood com reparametrizacao para ordenar alphas
neg_log_likelihood_ls <- function(params, data) {
  # params = (eta, a, delta, beta)
  eta   <- params[1]
  a     <- params[2]
  delta <- params[3]
  beta  <- params[4]

  pi1 <- exp(eta) / (1 + exp(eta))
  pi2 <- 1 - pi1

  alpha1 <- a
  alpha2 <- a + exp(delta) # garante alpha2 > alpha1

  f1 <- dgamma(data, shape = alpha1, rate = beta)
  f2 <- dgamma(data, shape = alpha2, rate = beta)
  density_mixture <- pi1 * f1 + pi2 * f2

  density_mixture[density_mixture < 1e-12] <- 1e-12
  -sum(log(density_mixture))
}
