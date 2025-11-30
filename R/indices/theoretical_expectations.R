# R/indices/theoretical_expectations.R
# Expectativas analíticas (E[estimadores]) sob a mistura Gamma.
# Depende de math_utils::compositions, log_mult_weights, softmax_from_log

compute_expected_value_TT_estimator <- function(n, pis, alphas) {
  m <- length(pis)
  stopifnot(length(alphas) == m, n >= 2)

  comps <- compositions(n, m)
  if (!is.matrix(comps)) comps <- t(as.matrix(comps))
  lgamma_n1 <- lgamma(n + 1)
  log_pis <- log(pis)

  log_p_k <- apply(comps, 1, function(k) lgamma_n1 - sum(lgamma(k + 1)) + sum(k * log_pis))
  log_p_k <- log_p_k - max(log_p_k)
  p_k <- exp(log_p_k); p_k <- p_k / sum(p_k)

  S_k <- as.vector(comps %*% alphas)
  trigamma_terms <- apply(comps, 1, function(k) sum(k * alphas * digamma(alphas)))
  digamma_terms <- digamma(S_k)
  term <- (trigamma_terms - S_k * digamma_terms + n - 1) / S_k
  sum(p_k * term) + log(n)
}

compute_expected_value_TL_estimator <- function(n, pis, alphas) {
  m <- length(pis)
  comps <- compositions(n, m)
  lgamma_n1 <- lgamma(n + 1)
  log_pis <- log(pis)
  ev_S <- 0
  for (r in seq_len(nrow(comps))) {
    ks <- comps[r, ]
    lp_k <- lgamma_n1 - sum(lgamma(ks + 1)) + sum(ks * log_pis)
    p_k <- exp(lp_k)
    S_k <- sum(ks * alphas)
    ev_S <- ev_S + p_k * digamma(S_k)
  }
  ev_S - sum(pis * digamma(alphas)) - log(n)
}

compute_expected_value_A1_estimator <- function(n, pis, alphas) {
  comps <- compositions(n, length(pis))
  lg_n1 <- lgamma(n + 1)
  log_pis <- log(pis)
  lg_ratio <- lgamma(alphas + 1 / n) - lgamma(alphas)
  ev_term <- 0
  for (r in seq_len(nrow(comps))) {
    K <- comps[r, ]
    lp_K <- lg_n1 - sum(lgamma(K + 1)) + sum(K * log_pis)
    S_K <- sum(K * alphas)
    log_dir_moment <- sum(K * lg_ratio)
    ev_term <- ev_term + exp(lp_K + log_dir_moment) / S_K
  }
  1 - n * ev_term
}

compute_expected_value_AInf_estimator <- function(n, pis, alphas,
                                                  rel.tol = 1e-9,
                                                  subdivisions = 3000L,
                                                  tail_eps = 1e-12) {
  # Implementação para m = 2 (como você já tinha)
  if (length(pis) != 2 || length(alphas) != 2) stop("AInf implementado apenas para mistura com m = 2.")
  pi1 <- pis[1]; pi2 <- pis[2]
  a1 <- alphas[1]; a2 <- alphas[2]
  lQ <- function(u, a) pgamma(u, shape = a, rate = 1, lower.tail = FALSE, log.p = TRUE)
  q_hi <- function(a) qgamma(1 - tail_eps, shape = a, rate = 1)
  Umax <- max(q_hi(a1), q_hi(a2)); if (!is.finite(Umax) || Umax <= 0) Umax <- 1e3
  b1 <- min(1, Umax / 10); b2 <- min(10, Umax / 2); segs <- unique(sort(c(0, b1, b2, Umax)))
  I_k <- function(k) {
    f <- function(u) exp(k * lQ(u, a1) + (n - k) * lQ(u, a2))
    val <- 0
    for (i in seq_len(length(segs) - 1L)) {
      part <- integrate(f, lower = segs[i], upper = segs[i + 1], rel.tol = rel.tol,
                        subdivisions = subdivisions, stop.on.error = FALSE)
      if (!is.null(attr(part, "message"))) {
        mid <- (segs[i] + segs[i + 1]) / 2
        p1 <- integrate(f, segs[i], mid, rel.tol = rel.tol, subdivisions = subdivisions, stop.on.error = FALSE)
        p2 <- integrate(f, mid, segs[i + 1], rel.tol = rel.tol, subdivisions = subdivisions, stop.on.error = FALSE)
        val <- val + p1$value + p2$value
      } else {
        val <- val + part$value
      }
    }
    val
  }
  k <- 0:n
  Sk <- a1 * k + a2 * (n - k)
  Ik <- vapply(k, I_k, numeric(1))
  w <- choose(n, k) * (pi1 ^ k) * (pi2 ^ (n - k))
  1 - n * sum(w * (Ik / Sk))
}

compute_expected_value_VMR_estimator <- function(n, pis, alphas, beta) {
  comps <- compositions(n, length(pis))
  if (!is.matrix(comps)) comps <- t(as.matrix(comps))
  lgamma_n1 <- lgamma(n + 1)
  log_pis <- log(pis)
  log_p_k <- apply(comps, 1, function(k) lgamma_n1 - sum(lgamma(k + 1)) + sum(k * log_pis))
  log_p_k <- log_p_k - max(log_p_k)
  p_k <- exp(log_p_k); p_k <- p_k / sum(p_k)
  S_k <- as.vector(comps %*% alphas)
  term1 <- 1 / (S_k + 1)
  term2 <- comps %*% (alphas * (alphas + 1))
  term3 <- S_k / n
  inner <- term1 * term2 - term3
  (n / ((n - 1) * beta)) * sum(p_k * inner)
}
