# R/indices/theoretical_indices.R
# Valores populacionais exatos para mistura Gamma
# Nota: este arquivo usa hypergeo::hypergeo quando necessário

# carregue hypergeo (já que não é um pacote aqui)
if (!requireNamespace("hypergeo", quietly = TRUE)) {
  stop("Instale o pacote 'hypergeo' para usar compute_Gini_G().")
}

library(hypergeo)

compute_Gini_G <- function(pii, alpha) {
  
  m <- length(pii)  #### Number of components
  
  #### Compute π*
  pi_star <- pii * alpha / sum(pii * alpha)
  
  #### Compute Beta function values
  B_ij <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      B_ij[i, j] <- beta(alpha[i], alpha[j] + 1)
    }
  }
  
  #### Compute 2F1 hypergeometric function at 1/2
  hypergeom_values <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      # ------------------------------------------------------------
      # Compute 2F1 hypergeometric function at 1/2
      # ------------------------------------------------------------
      # Nota: A função hypergeo() do pacote 'hypergeo' às vezes
      # retorna valores complexos mesmo quando o resultado esperado
      # é real (isso ocorre quando alguns parâmetros são negativos,
      # gerando pequenas partes imaginárias numéricas sem significado
      # estatístico). 
      # Aqui usamos as.numeric() para extrair apenas a parte real,
      # que é o valor correto do Gini teórico.
      hypergeom_values[i, j] <- as.numeric(hypergeo(alpha[i], -alpha[j], alpha[i] + 1, 1/2))
    }
  }
  
  #### Compute Gini coefficient G
  G <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      G <- G + (pii[i] * pi_star[j] * hypergeom_values[i, j] / (2^(alpha[i] - 1) * alpha[i] * B_ij[i, j]) )
    }
  }
  
  G <- G - 1  #### Apply the final adjustment in Equation (7)
  
  return(G)
}

compute_TheilT <- function(pii, alpha) {
  s <- sum(pii * alpha)
  sum_pi_alpha_psi <- sum(pii * alpha * digamma(alpha))
  (1 / s) * (sum_pi_alpha_psi + 1) - log(s)
}

compute_TheilL <- function(pii, alpha) {
  s <- sum(pii * alpha)
  sum_pi_psi <- sum(pii * digamma(alpha))
  log(s) - sum_pi_psi
}

compute_Atkinson1 <- function(pii, alpha) {
  s <- sum(pii * alpha)
  sum_pi_psi <- sum(pii * digamma(alpha))
  1 - (1 / s) * exp(sum_pi_psi)
}

compute_AtkinsonInf <- function(pii, alpha) {
  1
}

compute_VMR <- function(pii, alpha, beta) {
  alphapi <- sum(pii * alpha)
  alphapi_alpha1 <- sum(pii * alpha * (alpha + 1))
  (1 / (beta * alphapi)) * (alphapi_alpha1 - alphapi ^ 2)
}
