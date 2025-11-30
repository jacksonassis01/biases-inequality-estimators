# R/mixture/mixture_rng.R
# Gerador de amostras da mistura Gamma (2 componentes)

rMGamma <- function(n, pii, alpha, beta = 1) {
  if (length(pii) != length(alpha)) stop("pii e alpha devem ter mesmo comprimento")
  comps <- sample(seq_along(pii), size = n, replace = TRUE, prob = pii)
  y <- numeric(n)
  for (j in seq_along(pii)) {
    idx <- which(comps == j)
    if (length(idx) > 0) {
      y[idx] <- rgamma(length(idx), shape = alpha[j], rate = beta)
    }
  }
  y
}
