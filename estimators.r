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

Theil_T_estimator = function(X) {
  mu = mean(X)
  return(mean((X / mu) * log(X / mu)))
}

Theil_L_estimator = function(X) {
  mu = mean(X)
  return(mean(log(mu / X)))
}

Atkinson_1_estimator = function(X) {
  mu = mean(X)
  gm = exp(mean(log(X)))
  return(1 - gm / mu)
}

Atkinson_inf_estimator = function(X) {
  mu = mean(X)
  return(1 - min(X) / mu)
}

VMR_estimator = function(X) {
  mu = mean(X)
  return(var(X) / mu)
}
