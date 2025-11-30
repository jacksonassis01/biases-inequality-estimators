# main.R
# script principal para carregar o projeto reorganizado rapidamente

rm(list = ls())

# carrega todos os arquivos r
source("R/indices/math_utils.R")
source("R/indices/estimators.R")
source("R/indices/theoretical_indices.R")
source("R/indices/theoretical_expectations.R")
source("R/indices/biases.R")

source("R/mixture/mixture_loglik.R")
source("R/mixture/mixture_mle.R")
source("R/mixture/mixture_rng.R")

source("R/montecarlo/montecarlo_setup.R")
source("R/montecarlo/montecarlo_run.R")

# roda simulacao
set.seed(123)

results <- run_full_experiment(
  pii = c(0.6, 0.4),
  alpha1_values = c(0.5),
  alpha2_values = c(0.5, 1, 2, 3, 4, 5),
  beta = 1,
  n_values = c(seq(10, 50, 5), seq(60, 100, 10)),
  Nrep = 1000,
  verbose = TRUE
)

# calcula metricas mcse
results <- append_metrics(results)

# salva as simulacoes
# saveRDS(results, file = "results/R_1000_n_10_100_alpha1_0.50_alpha2_0.50_5_beta_1.00.rds")

# geracao das figuras e tabela mcse em dissertacao/
# aplicacao em dados reais em application-real-data/