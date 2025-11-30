# R/indices/math_utils.R
# Funções utilitárias usadas por índices e expectativas

compositions <- function(n, m) {
  if (m == 1) return(matrix(n, nrow = 1))
  do.call(rbind, lapply(0:n, function(i) cbind(i, compositions(n - i, m - 1))))
}

softmax_from_log <- function(logw) {
  lw <- logw - max(logw)
  w <- exp(lw)
  w / sum(w)
}

log_mult_weights <- function(K, pis, N) {
  # K: matrix (linhas = composições)
  lfactorial(N) - rowSums(lfactorial(K)) + (K %*% log(pis))[, 1]
}
