# R/montecarlo/montecarlo_setup.R
# Helpers e configuração para rodar Monte Carlo

# rpop: general wrapper para chamar geradores com lista de argumentos
rpop <- function(fpop, params) do.call(fpop, params)

# compute_indices: aplica lista de índices (cada elemento é lista com est_func, bias_func, true_func)
compute_indices <- function(sample, indices_list, est_params) {
  lapply(names(indices_list), function(name) {
    idx <- indices_list[[name]]
    theta_hat <- idx$est_func(sample)
    bias_est <- do.call(idx$bias_func, est_params)
    theta_corr <- theta_hat - bias_est
    list(
      index = name,
      theta_hat = theta_hat,
      theta_corr = theta_corr,
      bias_est = bias_est
    )
  })
}

# Exemplo de indices_list padrão (você pode sobrescrever ao chamar run_montecarlo)
# Cada entrada: list(est_func = function(x) ..., bias_func = function(n, pii, alpha, beta) ..., true_func = function(pii, alpha, beta) ...)
default_indices_list <- list(
  Gini = list(
    est_func = function(x) G_estimator(x),
    bias_func = function(n, pii, alpha, ...) compute_bias_Ghat_G_optimized(n, pii, alpha),
    true_func = function(pii, alpha, beta) compute_Gini_G(pii, alpha)
  ),
  TheilT = list(
    est_func = function(x) T_estimator(x, eps = 1),
    bias_func = function(n, pii, alpha, ...) compute_bias_TThat_TT(n, pii, alpha),
    true_func = function(pii, alpha, beta) compute_TheilT(pii, alpha)
  ),
  TheilL = list(
    est_func = function(x) T_estimator(x, eps = 0),
    bias_func = function(n, pii, alpha, ...) compute_bias_TLhat_TL(n, pii, alpha),
    true_func = function(pii, alpha, beta) compute_TheilL(pii, alpha)
  ),
  Atkinson1 = list(
    est_func = function(x) A_estimator(x, eps = 1),
    bias_func = function(n, pii, alpha, ...) compute_bias_A1hat_A1(n, pii, alpha),
    true_func = function(pii, alpha, beta) compute_Atkinson1(pii, alpha)
  ),
  AtkinsonInf = list(
    est_func = function(x) A_estimator(x, eps = Inf),
    bias_func = function(n, pii, alpha, ...) compute_bias_AInfhat_AInf(n, pii, alpha),
    true_func = function(pii, alpha, beta) compute_AtkinsonInf(pii, alpha)
  ),
  VMR = list(
    est_func = function(x) VMR_estimator(x),
    bias_func = function(n, pii, alpha, beta) compute_bias_VMRhat_VMR(n, pii, alpha, beta),
    true_func = function(pii, alpha, beta) compute_VMR(pii, alpha, beta)
  )
)
