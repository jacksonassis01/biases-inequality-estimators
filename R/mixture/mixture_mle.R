# R/mixture/mixture_mle.R
# Funções de estimação via MLE (usa neg_log_likelihood e neg_log_likelihood_ls)

# Assumimos que mixture_loglik.R já foi source() ou que está no mesmo ambiente.

estimate_gamma_mixture <- function(data, start_params = c(0.5, 2, 5, 1)) {
  # start_params: (pi1, alpha1, alpha2, beta)
  par0 <- c(log(start_params[1] / (1 - start_params[1])), start_params[2:4])

  mle_fit <- optim(
    par = par0,
    fn = neg_log_likelihood,
    data = data,
    method = "L-BFGS-B",
    lower = c(-Inf, 0.001, 0.001, 0.001),
    upper = c(Inf, Inf, Inf, Inf)
  )

  pi1_est <- exp(mle_fit$par[1]) / (1 + exp(mle_fit$par[1]))

  list(
    pi1 = pi1_est,
    pi2 = 1 - pi1_est,
    alpha1 = mle_fit$par[2],
    alpha2 = mle_fit$par[3],
    beta = mle_fit$par[4],
    logLik = -mle_fit$value,
    par_raw = mle_fit$par,
    optim = mle_fit
  )
}

estimate_gamma_mixture_ls <- function(data, start_params = c(0.5, 2, 1, 1)) {
  # start_params: (pi1, a, delta, beta)
  eta0 <- log(start_params[1] / (1 - start_params[1]))
  a0 <- start_params[2]
  delta0 <- start_params[3]
  beta0 <- start_params[4]

  par0 <- c(eta0, a0, delta0, beta0)

  mle_fit <- optim(
    par = par0,
    fn = neg_log_likelihood_ls,
    data = data,
    method = "L-BFGS-B",
    lower = c(-Inf, 0.001, log(1e-8), 0.001),
    upper = c(Inf, Inf, log(50), Inf)
  )

  eta_est <- mle_fit$par[1]
  a_est <- mle_fit$par[2]
  delta_est <- mle_fit$par[3]
  beta_est <- mle_fit$par[4]

  pi1_est <- exp(eta_est) / (1 + exp(eta_est))
  pi2_est <- 1 - pi1_est

  list(
    pi1 = pi1_est,
    pi2 = pi2_est,
    alpha1 = a_est,
    alpha2 = a_est + exp(delta_est),
    beta = beta_est,
    logLik = -mle_fit$value,
    par_raw = mle_fit$par,
    optim = mle_fit
  )
}
