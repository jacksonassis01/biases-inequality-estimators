rm(list = ls())

results <- readRDS("results/R_1000_n_10_100_alpha1_0.50_alpha2_0.50_5_beta_1.00.rds")

# indices selecionados
selected_indices <- c("TheilT", "TheilL", "Atkinson1", "AtkinsonInf", "VMR")

# figuras ======================================================================
source("dissertation/figures.R")

# grafico variando alpha2, fixando n = 15 e alpha1 = 0.5
fig_alpha2 <- plot_indices_by_alpha2_from_results(
  results = results,
  n_selected = 15,
  alpha1_selected = 0.5,
  indices = selected_indices
)

print(fig_alpha2)

# ggsave(
#   filename = "dissertation/figures/fig_vies_eqm_R_1000_a1_0.5_a2_0.5_5_n_15.pdf",
#   plot = fig_alpha2,
#   width = 10,
#   height = 14,
#   dpi = 300,
#   device = cairo_pdf
# )

# grafico variando n, fixando alpha1 = 0.5 e alpha2 = 2
fig_n <- plot_indices_by_n_from_results(
  results = results,
  alpha1_selected = 0.5,
  alpha2_selected = 2,
  indices = selected_indices
)

print(fig_n)

# ggsave(
#   filename = "dissertation/figures/fig_vies_eqm_R_1000_a1_0.5_a2_2_n_10_100.pdf",
#   plot = fig_n,
#   width = 10,
#   height = 14,
#   dpi = 300,
#   device = cairo_pdf
# )

#===============================================================================

# tabela mcse ==================================================================
source("dissertation/mcse_table.R")

# build_mcse_table(results, indices = selected_indices)

#===============================================================================