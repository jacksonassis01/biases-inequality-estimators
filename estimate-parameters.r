rm(list = ls())

set.seed(42)

x = rgamma(1000, 2, 1)
y = rgamma(1000, 3, 1)

z = c(x, y)

par(mfrow = c(1,2))

hist(z)

sample_mixed_gammas <- function(n, p = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  is_x1 <- runif(n) < p
  # gerar apenas o necessário:
  k <- sum(is_x1)
  samples <- numeric(n)
  samples[is_x1] <- rgamma(k, shape = 2, scale = 1)
  samples[!is_x1] <- rgamma(n - k, shape = 3, scale = 1)
  return(samples)
}

z2 <- sample_mixed_gammas(2000, p = 0.5, seed = 42)

hist(z2)

fit = gammamixEM(z)

fit$lambda
fit$gamma.pars

fit2 = gammamixEM(z2)

fit2$lambda
fit2$gamma.pars
