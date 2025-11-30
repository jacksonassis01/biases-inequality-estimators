# R/indices/biases.R
# Viés analítico: bias = E[est] - true_index
# Depende de theoretical_expectations.R e theoretical_indices.R (ambos em indices/)
# Depende de math_utils.R em indices/

compute_bias_TThat_TT <- function(n, pii, alpha, ...) {
  compute_expected_value_TT_estimator(n, pii, alpha) - compute_TheilT(pii, alpha)
}

compute_bias_TLhat_TL <- function(n, pii, alpha, ...) {
  compute_expected_value_TL_estimator(n, pii, alpha) - compute_TheilL(pii, alpha)
}

compute_bias_A1hat_A1 <- function(n, pii, alpha, ...) {
  compute_expected_value_A1_estimator(n, pii, alpha) - compute_Atkinson1(pii, alpha)
}

compute_bias_AInfhat_AInf <- function(n, pii, alpha, ...) {
  compute_expected_value_AInf_estimator(n, pii, alpha) - compute_AtkinsonInf(pii, alpha)
}

compute_bias_VMRhat_VMR <- function(n, pii, alpha, beta) {
  compute_expected_value_VMR_estimator(n, pii, alpha, beta) - compute_VMR(pii, alpha, beta)
}

# Versão otimizada do viés do Gini (vetorizada, usa math_utils)
compute_bias_Ghat_G_optimized <- function(n, pii, alpha) {
  stopifnot(length(pii) == length(alpha), n >= 2)
  m <- length(pii)
  alpha_bar <- sum(pii * alpha)
  pi_star <- pii * alpha / alpha_bar

  comps_nm2 <- compositions(n - 2, m)
  if (!is.matrix(comps_nm2)) comps_nm2 <- t(as.matrix(comps_nm2))
  comps_nm1 <- compositions(n - 1, m)
  if (!is.matrix(comps_nm1)) comps_nm1 <- t(as.matrix(comps_nm1))

  log_p_k_2 <- log_mult_weights(comps_nm2, pii, n - 2)
  p_k_2 <- softmax_from_log(log_p_k_2)
  S_k_2 <- as.vector(comps_nm2 %*% alpha)

  log_p_k_1 <- log_mult_weights(comps_nm1, pii, n - 1)
  p_k_1 <- softmax_from_log(log_p_k_1)
  S_k_1 <- as.vector(comps_nm1 %*% alpha)

  coeff_ij <- function(i, j) 2 * pii[i] * pi_star[j] * pbeta(0.5, alpha[i], alpha[j] + 1)

  term1_total <- 0
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      denom_vec <- S_k_2 + alpha[i] + alpha[j]
      num_exp <- sum(p_k_2 * (n * alpha_bar) / denom_vec)
      adj <- num_exp - 1
      term1_total <- term1_total + coeff_ij(i, j) * adj
    }
  }

  term2_mat <- sweep(matrix(S_k_1, nrow = length(S_k_1), ncol = m), 2, alpha, `+`)
  numer_j <- matrix(rep(n * (pii * alpha), each = nrow(term2_mat)), nrow = nrow(term2_mat), ncol = m)
  term2_inner <- rowSums(numer_j / term2_mat)
  Bias2_G <- 1 - sum(p_k_1 * term2_inner)

  term1_total + Bias2_G
}

# Correção aplicada ao estimador amostral
bias_corrected_Ghat <- function(x, pii, alpha) {
  G_estimator(x) - compute_bias_Ghat_G_optimized(length(x), pii, alpha)
}

bias_corrected_That <- function(x, pii, alpha) {
  T_estimator(x, eps = 1) - compute_bias_TThat_TT(length(x), pii, alpha)
}

bias_corrected_Lhat <- function(x, pii, alpha) {
  T_estimator(x, eps = 0) - compute_bias_TLhat_TL(length(x), pii, alpha)
}

bias_corrected_A1hat <- function(x, pii, alpha) {
  A_estimator(x, eps = 1) - compute_bias_A1hat_A1(length(x), pii, alpha)
}

bias_corrected_AInfhat <- function(x, pii, alpha) {
  A_estimator(x, eps = Inf) - compute_bias_AInfhat_AInf(length(x), pii, alpha)
}

bias_corrected_VMRhat <- function(x, pii, alpha, beta) {
  VMR_estimator(x) - compute_bias_VMRhat_VMR(length(x), pii, alpha, beta)
}
