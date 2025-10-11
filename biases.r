library(gtools)
library(rlist)
library(hypergeo)

compositions = function(n, m) {
  if (m == 1)
    return(matrix(n, nrow = 1))
  do.call(rbind, lapply(0:n, function(i)
    cbind(i, compositions(n - i, m - 1))))
}

# compositions = function(n, m) {
#   if (m == 1)
#     return(matrix(n, nrow = 1))
#   res = NULL
#   for (i in 0:n) {
#     sub = compositions(n - i, m - 1)
#     res = rbind(res, cbind(i, sub))
#   }
#
#   res
# }

# gini =========================================================================
compute_bias_Ghat_G <- function(n, pii, alpha, ...) {
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
      hypergeom_values[i, j] <- as.numeric(hypergeo(-alpha[j], alpha[i], alpha[i] + 1, 1 / 2))
    }
  }
  
  #### Compute First Bias Term
  Bias1_G_hat_G <- 0
  
  for (i in 1:m) {
    for (j in 1:m) {
      term1 <- 0
      
      #### Generate all index combinations for (n-2) elements
      subset_indices <- expand.grid(rep(list(1:m), (n - 2)))
      
      for (k in 1:nrow(subset_indices)) {
        subset_vector <- as.numeric(subset_indices[k, ])  #### Convert to numeric vector
        subset_alpha  <- sum(alpha[subset_vector])
        subset_pi_prod <- prod(pii[subset_vector])
        term1 <- term1 + ((subset_pi_prod * n * sum(pii * alpha)) / (subset_alpha + alpha[i] + alpha[j]))
      }
      
      term1 <- term1 - 1
      
      Bias1_G_hat_G <- Bias1_G_hat_G + (pii[i] * pi_star[j] * hypergeom_values[i, j] /
                                          (2 ^ (alpha[i] - 1) * alpha[i] * B_ij[i, j])) * term1
    }
  }
  
  #### Compute Second Bias Term
  Bias2_G_hat_G <- 0
  
  #### Generate all index combinations for (n-1) elements
  subset_indices <- expand.grid(rep(list(1:m), (n - 1)))
  
  for (k in 1:nrow(subset_indices)) {
    subset_vector <- as.numeric(subset_indices[k, ])  #### Convert to numeric vector
    subset_alpha <- sum(alpha[subset_vector])
    subset_pi_prod <- prod(pii[subset_vector])
    term2 <- 0
    for (j in 1:m) {
      term2 <- term2 + (n * pii[j] * alpha[j]) / (subset_alpha + alpha[j])
    }
    Bias2_G_hat_G <- Bias2_G_hat_G + subset_pi_prod * term2
  }
  
  Bias2_G_hat_G <- 1 - Bias2_G_hat_G
  
  #### Final Bias Value
  Bias_G_hat_G <- Bias1_G_hat_G + Bias2_G_hat_G
  
  return(Bias_G_hat_G)
}

compute_bias_Ghat_G_fast <- function(n, pii, alpha, ...) {
  if (!requireNamespace("hypergeo", quietly = TRUE))
    stop("O pacote 'hypergeo' é necessário. Instale-o com install.packages('hypergeo').")
  
  m <- length(pii)
  alpha_bar <- sum(pii * alpha)
  pi_star   <- pii * alpha / alpha_bar
  
  # Matrizes auxiliares (i,j)
  B_ij <- outer(alpha, alpha, function(a_i, a_j)
    beta(a_i, a_j + 1))
  H_ij <- outer(alpha, alpha, Vectorize(function(a_i, a_j) {
    as.numeric(hypergeo::hypergeo(-a_j, a_i, a_i + 1, 0.5))
  }))
  
  # -------- Termo 1: somas sobre (n-2) --------
  comps_nm2 <- compositions(n - 2, m)
  fact_nm2  <- factorial(n - 2)
  fact_ks2  <- apply(comps_nm2, 1, function(k)
    prod(factorial(k)))
  pii_pow2  <- t(apply(comps_nm2, 1, function(k)
    pii ^ k))
  p_k_2     <- fact_nm2 / fact_ks2 * apply(pii_pow2, 1, prod)   # pesos multinomiais de (n-2)
  S_k_2     <- as.vector(comps_nm2 %*% alpha)                    # sum_j k_j * alpha_j (n-2)
  
  term1_total <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      # média ponderada correta sobre (n-2)
      num <- sum(p_k_2 * (n * alpha_bar) / (S_k_2 + alpha[i] + alpha[j]))
      adj <- num - 1
      coeff <- (pii[i] * pi_star[j] * H_ij[i, j]) /
        (2 ^ (alpha[i] - 1) * alpha[i] * B_ij[i, j])
      term1_total <- term1_total + coeff * adj
    }
  }
  
  # -------- Termo 2: somas sobre (n-1) --------
  comps_nm1 <- compositions(n - 1, m)
  fact_nm1  <- factorial(n - 1)
  fact_ks1  <- apply(comps_nm1, 1, function(k)
    prod(factorial(k)))
  pii_pow1  <- t(apply(comps_nm1, 1, function(k)
    pii ^ k))
  p_k_1     <- fact_nm1 / fact_ks1 * apply(pii_pow1, 1, prod)   # pesos multinomiais de (n-1)
  S_k_1     <- as.vector(comps_nm1 %*% alpha)                    # sum_j k_j * alpha_j (n-1)
  
  term2_inner <- sapply(seq_along(S_k_1), function(idx) {
    sum(n * pii * alpha / (S_k_1[idx] + alpha))
  })
  Bias2_G <- 1 - sum(p_k_1 * term2_inner)
  
  # -------- Viés total --------
  Bias_Ghat_G <- term1_total + Bias2_G
  return(Bias_Ghat_G)
}

#===============================================================================

# theil t ======================================================================
compute_bias_theil_t_hat_theil_t_1 = function(n, pii, alpha, ...) {
  m = length(pii)
  
  indices = permutations(
    n = m,
    r = n,
    v = 1:m,
    repeats.allowed = TRUE
  )
  
  sum1 = 0
  for (k in 1:nrow(indices)) {
    j = indices[k, ]
    pi_prod = prod(pii[j])
    alpha_sum = sum(alpha[j])
    inside = sum(alpha[j] * digamma(alpha[j])) - alpha_sum * digamma(alpha_sum) + n - 1
    sum1 = sum1 + pi_prod * (inside / alpha_sum)
  }
  
  sum1 = sum1 + log(n)
  
  alpha_bar = sum(pii * alpha)
  sum2 = (sum(pii * alpha * digamma(alpha)) + 1) / alpha_bar - log(alpha_bar)
  
  bias = sum1 - sum2
  
  return(bias)
}

compute_bias_theil_t_hat_theil_t = function(n, pii, alpha, ...) {
  comps = expand.grid(rep(list(1:length(pii)), n))
  
  E_T_T_hat = 0
  for (i in 1:nrow(comps)) {
    js = as.numeric(comps[i, ])
    k = prod(pii[js])
    s = sum(alpha[js])
    term = sum(alpha[js] * digamma(alpha[js])) - (s * digamma(s)) + n -
      1
    E_T_T_hat = E_T_T_hat + (k * (1 / s) * term)
  }
  
  E_T_T_hat = E_T_T_hat + log(n)
  
  T_T = (1 / sum(pii * alpha)) * (sum(pii * alpha * digamma(alpha)) + 1) - log(sum(pii *
                                                                                     alpha))
  
  bias = E_T_T_hat - T_T
  
  return(bias)
}

compute_bias_TThat_TT_fast = function(n, pii, alpha, ...) {
  m = length(pii)
  
  comps = compositions(n, m)
  fact_n = factorial(n)
  
  fact_ks = apply(comps, 1, function(k)
    prod(factorial(k)))
  pii_pow = t(apply(comps, 1, function(k)
    pii ^ k))
  p_k = fact_n / fact_ks * apply(pii_pow, 1, prod)
  
  S_k = comps %*% alpha
  
  trigamma_terms = apply(comps, 1, function(k)
    sum(k * alpha * digamma(alpha)))
  digamma_terms = digamma(S_k)
  term = (trigamma_terms - S_k * digamma_terms + n - 1) / S_k
  
  E_T_T_hat = sum(p_k * term) + log(n)
  
  T_T = (1 / sum(pii * alpha)) * (sum(pii * alpha * digamma(alpha)) + 1) - log(sum(pii * alpha))
  
  bias = E_T_T_hat - T_T
  
  return(bias)
}

n = 15
pii = c(.45, .55)
alpha = c(.5, 3)
compute_bias_theil_t_hat_theil_t_1(n, pii, alpha)
compute_bias_theil_t_hat_theil_t(n, pii, alpha)
compute_bias_TThat_TT_fast(n, pii, alpha)

#===============================================================================

# theil l ======================================================================
compute_bias_theil_l_hat_theil_l_1 = function(n, pii, alpha, ...) {
  m = length(pii)
  if (length(alpha) != m)
    stop("pii and alpha must have same length")
  if (abs(sum(pii) - 1) > 1e-12)
    stop("pii must sum to 1")
  
  alpha_bar = sum(pii * alpha)            # sum_j pi_j * alpha_j
  
  # gera todas as composições k1..km que somam n
  compositions = function(n, m) {
    if (m == 1)
      return(matrix(n, nrow = 1))
    res = NULL
    for (i in 0:n) {
      sub = compositions(n - i, m - 1)
      res = rbind(res, cbind(i, sub))
    }
    res
  }
  
  comps = compositions(n, m)   # cada linha = (k1,..,km)
  fact_n = factorial(n)
  big_sum = 0
  for (r in seq_len(nrow(comps))) {
    ks = comps[r, ]
    coeff = fact_n / prod(factorial(ks))          # coef. multinomial
    p_k = coeff * prod(pii ^ ks)                  # P(k1,..,km)
    S_k = sum(ks * alpha)
    big_sum = big_sum + p_k * digamma(S_k)        # E[ psi(S) ]
  }
  
  bias_TL = big_sum - log(n) - log(alpha_bar)
  return(bias_TL)
}

compute_bias_theil_l_hat_theil_l_2 = function(n, pii, alpha, ...) {
  m = length(pii)
  
  alpha_bar = sum(pii * alpha)            # sum_j pi_j * alpha_j
  
  # gera todas as composições k1..km que somam n
  compositions = function(n, m) {
    if (m == 1)
      return(matrix(n, nrow = 1))
    res = NULL
    for (i in 0:n) {
      sub = compositions(n - i, m - 1)
      res = rbind(res, cbind(i, sub))
    }
    res
  }
  
  comps = compositions(n, m)   # cada linha = (k1,..,km)
  fact_n = factorial(n)
  big_sum = 0
  for (r in seq_len(nrow(comps))) {
    ks = comps[r, ]
    coeff = fact_n / prod(factorial(ks))          # coef. multinomial
    p_k = coeff * prod(pii ^ ks)                  # P(k1,..,km)
    S_k = sum(ks * alpha)
    big_sum = big_sum + p_k * digamma(S_k)        # E[ psi(S) ]
  }
  
  bias_TL = big_sum - log(n) - log(alpha_bar)
  return(bias_TL)
}

compute_bias_theil_l_hat_theil_l_3 = function(n, pii, alpha, ...) {
  comps = expand.grid(rep(list(1:length(pii)), n))
  
  ms = 0
  
  for (i in 1:nrow(comps)) {
    js = as.numeric(comps[i, ])
    ms = ms + (prod(pii[js]) * digamma(sum(alpha[js])))
  }
  
  bias = ms - log(n) - log(sum(pii * alpha))
  
  return(bias)
}

compute_bias_theil_l_hat_theil_l = function(n, pii, alpha, ...) {
  comps = expand.grid(rep(list(1:length(pii)), n))
  
  ms = 0
  
  for (i in 1:nrow(comps)) {
    js = as.numeric(comps[i, ])
    ms = ms + (prod(pii[js]) * digamma(sum(alpha[js])))
  }
  
  bias = ms - log(n) - log(sum(pii * alpha))
  
  return(bias)
}

compute_bias_TLhat_TL_fast = function(n, pii, alpha, ...) {
  m = length(pii)
  
  comps = compositions(n, m)
  
  fact_n = factorial(n)
  big_sum = 0
  for (r in seq_len(nrow(comps))) {
    ks = comps[r, ]
    coeff = fact_n / prod(factorial(ks))
    p_k = coeff * prod(pii ^ ks)
    S_k = sum(ks * alpha)
    big_sum = big_sum + p_k * digamma(S_k)
  }
  
  bias_TL = big_sum - log(n) - log(sum(pii * alpha))
  return(bias_TL)
}

n = 15
pii = c(.45, .55)
alpha = c(.5, 3)
compute_bias_theil_l_hat_theil_l_1(n, pii, alpha)
compute_bias_theil_l_hat_theil_l_2(n, pii, alpha)
compute_bias_theil_l_hat_theil_l_3(n, pii, alpha)
compute_bias_theil_l_hat_theil_l(n, pii, alpha)
compute_bias_TLhat_TL_fast(n, pii, alpha)

#===============================================================================

# atkinson 1 ===================================================================
compute_bias_A1hat_A1 = function(n, pii, alpha, ...) {
  if (abs(sum(pii) - 1) > 1e-12)
    stop("pii must sum to 1")
  m = length(pii)
  if (length(alpha) != m)
    stop("pii and alpha must have same length")
  alpha_bar = sum(pii * alpha)
  term2 = (1 / alpha_bar) * exp(sum(pii * digamma(alpha)))
  # gera composições k1..km que somam n
  compositions = function(n, m) {
    if (m == 1)
      return(matrix(n, nrow = 1))
    res = NULL
    for (i in 0:n) {
      sub = compositions(n - i, m - 1)
      res = rbind(res, cbind(i, sub))
    }
    res
  }
  comps = compositions(n, m)
  big_sum = 0
  fact_n = factorial(n)
  for (r in seq_len(nrow(comps))) {
    ks = comps[r, ]
    coeff = fact_n / prod(factorial(ks))
    p_k = coeff * prod(pii ^ ks)
    S_k = sum(ks * alpha)
    # produto (Gamma(alpha_j + 1/n) / Gamma(alpha_j))^k_j
    prod_factor = prod((gamma(alpha + 1.0 / n) / gamma(alpha)) ^ ks)
    big_sum = big_sum + p_k * (prod_factor / S_k)
  }
  bias = -n * big_sum + term2
  return(bias)
}

#===============================================================================

# atkinson inf =================================================================
compute_bias_Ainfhat_Ainf_1 = function(n,
                                       pii,
                                       alpha,
                                       u_max = 100,
                                       du = 0.1,
                                       ...) {
  # Gera composições fracas: todas as tuplas (k1, ..., km) com soma n
  generate_multinomial_indices <- function(n, m) {
    result <- list()
    helper <- function(prefix, remaining, depth) {
      if (depth == m) {
        result[[length(result) + 1]] <<- c(prefix, remaining)
        return(invisible(NULL))
      }
      for (i in 0:remaining) {
        helper(c(prefix, i), remaining - i, depth + 1)
      }
    }
    helper(c(), n, 1)
    do.call(rbind, result)
  }
  
  m = length(pii)
  
  # Gera todas as composições fracas de n em m partes
  partitions = generate_multinomial_indices(n, m)
  
  bias_sum = 0
  
  for (row in 1:nrow(partitions)) {
    k = partitions[row, ]
    
    # Coeficiente multinomial: n! / (k1! * k2! * ... * km!)
    multinom_coeff = factorial(n) / prod(factorial(k))
    
    # Produto das pi^k
    pi_term = prod(pii ^ k)
    
    # Soma ponderada dos alpha
    alpha_dot_k = sum(alpha * k)
    if (alpha_dot_k == 0)
      next  # evitar divisão por zero
    
    # Função integrando
    integrand = function(u) {
      prod_term = 1
      for (j in 1:m) {
        if (k[j] == 0) {
          term = 1
        } else {
          # Gama incompleta superior: Γ(a, u) = pgamma(u, a, lower.tail = FALSE) * gamma(a)
          gamma_inc = pgamma(u, shape = alpha[j], lower.tail = FALSE) * gamma(alpha[j])
          term = (gamma_inc / gamma(alpha[j])) ^ k[j]
        }
        prod_term = prod_term * term
      }
      return(prod_term)
    }
    
    # Integral numérica via soma de Riemann (pode substituir por integrate se preferir)
    u_vals = seq(0, u_max, by = du)
    integrand_vals = sapply(u_vals, integrand)
    integral_approx = sum(integrand_vals) * du
    
    term_k = multinom_coeff * pi_term * (1 / alpha_dot_k) * integral_approx
    
    bias_sum = bias_sum + term_k
  }
  
  bias = -n * bias_sum
  return(bias)
}

compute_bias_Ainfhat_Ainf = function(n, pii, alpha, ...) {
  q = function(a, u) {
    return(pgamma(u, shape = a, lower.tail = FALSE))
  }
  
  qq = function(a1, a2, n, i, u) {
    return((q(a1, u) ^ i / gamma(a1)) * (q(a2, u) ^ (n - i) / gamma(a2)))
  }
  
  s = 0
  
  for (i in 0:n) {
    f = factorial(n) / (factorial(i) * factorial(n - i))
    p = pii[1] ^ i * pii[2] ^ (n - i)
    inv = 1 / (alpha[1] * i + alpha[1] * (n - i))
    I = integrate(\(u) qq(alpha[1], alpha[2], n, i, u), 50, Inf)$value +
      integrate(\(u) qq(alpha[1], alpha[2], n, i, u), 50, Inf)$value
    s = s + f * p * inv * I
  }
  
  return(-n * s)
}

compute_bias_Ainfhat_Ainf_fast <- function(n,
                                           pii,
                                           alpha,
                                           u_max = 100,
                                           du = 0.1,
                                           ...) {
  m <- length(pii)
  K  <- compositions(n, m)                # cada linha: k = (k1,...,km)
  nK <- nrow(K)
  
  # coeficientes multinomiais e termos de prob
  fact_n  <- factorial(n)
  fact_ks <- apply(K, 1, function(k)
    prod(factorial(k)))
  w_k     <- fact_n / fact_ks * apply(K, 1, function(k)
    prod(pii ^ k))  # pesos p_k
  denom   <- as.vector(K %*% alpha)                                   # alpha·k
  keep    <- denom > 0
  if (!any(keep))
    return(0)
  
  # grade para integrar (pode ajustar u_max/du conforme precisão)
  u_vals <- seq(0, u_max, by = du)
  
  # matriz log Q(u; alpha_j): len(u) x m
  # cuidado com log(0): substitui -Inf por um valor muito negativo para não gerar NaN
  logQ <- sapply(alpha, function(a) {
    v <- pgamma(u_vals, shape = a, lower.tail = FALSE)
    lv <- log(v)
    lv[!is.finite(lv)] <- -1e300
    lv
  })
  
  # integral por composição: ∑_u exp( (logQ %*% k) ) * du
  integrals <- apply(K[keep, , drop = FALSE], 1, function(k) {
    sum(exp(logQ %*% k)) * du
  })
  
  terms <- w_k[keep] * (integrals / denom[keep])
  - n * sum(terms[is.finite(terms)])
}

n = 15
pii = c(.45, .55)
alpha = c(.5, 3)
compute_bias_Ainfhat_Ainf_1(n, pii, alpha)
compute_bias_Ainfhat_Ainf(n, pii, alpha)
compute_bias_Ainfhat_Ainf_fast(n, pii, alpha)

#===============================================================================

# vmr ==========================================================================
compute_bias_vmr_hat_vmr = function(n, pii, alpha, beta, ...) {
  m = length(pii)
  
  indices = permutations(
    n = m,
    r = n,
    v = 1:m,
    repeats.allowed = TRUE
  )
  
  sum1 = 0
  for (k in 1:nrow(indices)) {
    js = indices[k, ]
    pii_prod = prod(pii[js])
    
    alpha_js = alpha[js]
    alpha_sum = sum(alpha_js)
    
    term1 = 1 / (alpha_sum + 1)
    term2 = sum(alpha_js * (alpha_js + 1))
    term3 = alpha_sum / n
    
    inner_term = (term1 * term2) - term3
    
    sum1 = sum1 + (pii_prod * inner_term)
  }
  
  first_part = (n / ((n - 1) * beta)) * sum1
  
  alpha_bar = sum(pii * alpha)
  term_num = sum(pii * alpha * (alpha + 1)) - (alpha_bar) ^ 2
  second_part = (1 / (beta * alpha_bar)) * term_num
  
  bias = first_part - second_part
  
  return(bias)
}

compute_bias_VMRhat_VMR_fast = function(n, pii, alpha, beta, ...) {
  m = length(pii)
  
  comps = compositions(n, m)
  fact_n = factorial(n)
  
  fact_ks = apply(comps, 1, function(k)
    prod(factorial(k)))
  pii_pow = t(apply(comps, 1, function(k)
    pii ^ k))
  p_k = fact_n / fact_ks * apply(pii_pow, 1, prod)
  
  S_k = comps %*% alpha
  
  term1 = 1 / (S_k + 1)
  term2 = comps %*% (alpha * (alpha + 1))
  term3 = S_k / n
  inner = term1 * term2 - term3
  
  first_part = (n / ((n - 1) * beta)) * sum(p_k * inner)
  
  alpha_bar = sum(pii * alpha)
  num = sum(pii * alpha * (alpha + 1)) - (alpha_bar) ^ 2
  second_part = (1 / (beta * alpha_bar)) * num
  
  bias = first_part - second_part
  
  return(bias)
}

n = 15
pii = c(.45, .55)
alpha = c(.5, 3)
beta = 1
compute_bias_vmr_hat_vmr(n, pii, alpha, beta)
compute_bias_VMRhat_VMR_fast(n, pii, alpha, beta)
#===============================================================================