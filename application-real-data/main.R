#===============================================================================
# aplicacao empirica: ajuste da mistura gama em dados de pib continentais
#===============================================================================

rm(list = ls())

library(dplyr)
library(ggplot2)
library(scales)

source("R/mixture/mixture_loglik.R")
source("R/mixture/mixture_mle.R")
source("R/indices/estimators.R")
source("R/indices/math_utils.R")
source("R/indices/theoretical_indices.R")
source("R/indices/theoretical_expectations.R")
source("R/indices/biases.R")

# importacao dos dados e estatisticas descritivas ==============================
df <- read.csv("application-real-data/data/processed/gdp-per-capita-worldbank-with-continent.csv") %>%
  filter(Continent %in% c("North America", "Oceania"))

X <- df$gdp_per_capita
n <- nrow(df)

# medidas resumo
media <- mean(X)
mediana <- median(X)
desvio_padrao <- sd(X)

medidas_resumo <- data.frame(
  Medida = c("Media", "Mediana", "Desvio-Padrao"),
  Valor = c(media, mediana, desvio_padrao)
)

# arredondamento para 4 casas decimais na exibicao
medidas_resumo$Valor <- round(medidas_resumo$Valor, 4)

print("--- Medidas Resumo (mil dolares) ---")
print(medidas_resumo)

# dados para o histograma
hist_data <- hist(X, plot = FALSE)
df_hist <- data.frame(
  x = hist_data$mids,
  count = hist_data$counts,
  density = hist_data$density
)

#===============================================================================

# ajuste da mistura gama de 2 componentes ======================================
fit <-  estimate_gamma_mixture_ls(X)

# Parametros estimados
pii <- c(fit$pi1, fit$pi2)
alphas <- c(fit$alpha1, fit$alpha2)
beta <- fit$beta

# Parametros arredondados para exibicao
fit_params <- c(pii, alphas, beta)
names(fit_params) <- c("pi1", "pi2", "alpha1", "alpha2", "lambda")
fit_params_rounded <- round(fit_params, 4)

print("--- Parametros da Mistura Gama ---")
print(fit_params_rounded)

#===============================================================================

# calculo dos indices de desigualdade e correcoes de vies ======================
indices <- list(
  Gini = list(
    est_func = function(x) G_estimator(x),
    bias_func = function(n, pii, alpha, ...) compute_bias_Ghat_G_optimized(n, pii, alpha)
  ),
  TheilT = list(
    est_func = function(x) T_estimator(x, eps = 1),
    bias_func = function(n, pii, alpha, ...) compute_bias_TThat_TT(n, pii, alpha)
  ),
  TheilL = list(
    est_func = function(x) T_estimator(x, eps = 0),
    bias_func = function(n, pii, alpha, ...) compute_bias_TLhat_TL(n, pii, alpha)
  ),
  Atkinson1 = list(
    est_func = function(x) A_estimator(x, eps = 1),
    bias_func = function(n, pii, alpha, ...) compute_bias_A1hat_A1(n, pii, alpha)
  ),
  AtkinsonInf = list(
    est_func = function(x) A_estimator(x, eps = Inf),
    bias_func = function(n, pii, alpha, ...) compute_bias_AInfhat_AInf(n, pii, alpha)
  ),
  VMR = list(
    est_func = function(x) VMR_estimator(x),
    bias_func = function(n, pii, alpha, beta) compute_bias_VMRhat_VMR(n, pii, alpha, beta)
  )
)

results <- list()

for (idx in names(indices)) {
  f <- indices[[idx]]
  
  theta_hat <- f$est_func(X)
  bias_est <- f$bias_func(n, pii, alphas, beta)
  
  theta_corr <- theta_hat - bias_est
  
  results[[idx]] <- data.frame(
    indice = idx,
    theta_hat = theta_hat,
    theta_corr = theta_corr,
    bias_est = bias_est
  )
}

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL

# arredondamento para 4 casas decimais nas colunas numericas
results_df$theta_hat <- round(results_df$theta_hat, 4)
results_df$theta_corr <- round(results_df$theta_corr, 4)
results_df$bias_est <- round(results_df$bias_est, 4)

print("--- Resultados dos Indices de Desigualdade ---")
print(results_df)

#===============================================================================

# figuras ======================================================================

# histograma de frequencia
g1 <- ggplot(df_hist, aes(x = x, y = count)) +
  geom_col(
    fill = "white",
    color = "black",
    width = diff(hist_data$breaks)[1]
  ) +
  scale_x_continuous(labels = label_number(suffix = "k")) +
  scale_y_continuous(
    breaks = seq(0, max(df_hist$count), by = 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "PIB per capita (mil dolares)", y = "Frequencia de paises") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_line(color = "gray96"),
    panel.grid.minor = element_line(color = "gray96")
  )

print(g1)

# ggsave(
#   filename = "application-real-data/figures/hist_pib_global_north_america_oceania.pdf",
#   plot = g1,
#   width = 10,
#   height = 14,
#   dpi = 300,
#   device = cairo_pdf
# )

# histograma de densidade com curva ajustada
dMixGamma <- function(x, pi1, pi2, alpha1, alpha2, beta) {
  pi1 * dgamma(x, shape = alpha1, rate = beta) +
    pi2 * dgamma(x, shape = alpha2, rate = beta)
}

x_grid <- seq(0, max(X) * 1.1, length.out = 1000)
y_mix <- dMixGamma(x_grid, pii[1], pii[2], alphas[1], alphas[2], beta)

df_curve <- data.frame(x = x_grid, y = y_mix)

g2 <- ggplot(df_hist, aes(x = x, y = density)) +
  geom_col(
    fill = NA,
    color = "black",
    linetype = "dashed",
    width = diff(hist_data$breaks)[1]
  ) +
  geom_line(
    data = df_curve,
    aes(x = x, y = y),
    color = "black",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 0) +
  scale_x_continuous(labels = label_number(suffix = "k")) +
  labs(x = "PIB per capita (mil dolares)", y = "Densidade") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_line(color = "gray96"),
    panel.grid.minor = element_line(color = "gray96")
  )

print(g2)

# ggsave(
#   filename = "application-real-data/figures/hist_pib_global_north_america_oceania_with_curve.pdf",
#   plot = g2,
#   width = 10,
#   height = 14,
#   dpi = 300,
#   device = cairo_pdf
# )

#===============================================================================
