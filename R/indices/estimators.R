# R/indices/estimators.R
# Estimadores amostrais (mantive nomes originais)

G_estimator <- function(x) {
  n <- length(x)
  if (n < 2) stop("Sample size must be at least 2.")
  num <- sum(abs(outer(x, x, "-"))) / 2
  den <- sum(x)
  (1 / (n - 1)) * (num / den)
}

T_estimator <- function(x, eps = 1) {
  if (!(eps %in% c(1, 0))) stop("eps deve ser 1 (Theil-T) ou 0 (Theil-L).")
  mu <- mean(x)
  if (eps == 1) {
    return(mean((x / mu) * log(x / mu)))
  } else {
    return(mean(log(mu / x)))
  }
}

A_estimator <- function(x, eps = 1) {
  if (!(eps %in% c(1, Inf))) stop("eps deve ser 1 ou Inf.")
  mu <- mean(x)
  if (eps == 1) {
    gm <- exp(mean(log(x)))
    return(1 - gm / mu)
  } else {
    return(1 - min(x) / mu)
  }
}

VMR_estimator <- function(x) {
  var(x) / mean(x)
}
