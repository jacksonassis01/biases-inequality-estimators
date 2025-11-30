#===============================================================================
# mcse_table.R
# Construcao de tabela latex para os Erros Padrao de Monte Carlo
#===============================================================================

library(dplyr)
library(stringr)

# funcao para converter o cenario em linhas latex ==============================

converte_cenario_em_latex <- function(lines) {
  lines_latex = ""
  for (i in 1:nrow(lines)) {
    if (i == 1) {
      s = sprintf(
        "\\multirow{24}{*}{%s} & \\multirow{4}{*}{%.1f} & \\multirow{4}{*}{%.4f} & %d & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\ \n",
        lines[i, 1], lines[i, 2], lines[i, 3], lines[i, 4], lines[i, 5], lines[i, 6],
        lines[i, 7], lines[i, 8], lines[i, 9], lines[i, 10], lines[i, 11], lines[i, 12],
        lines[i, 13], lines[i, 14]
      )
    } else if (i %in% (4 * (1:6) + 1)) {
      s = sprintf(
        "\n \\cmidrule(lr){2-14} & \\multirow{4}{*}{%.1f} & \\multirow{4}{*}{%.4f} & %d & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\ \n",
        lines[i, 2], lines[i, 3], lines[i, 4], lines[i, 5], lines[i, 6], lines[i, 7],
        lines[i, 8], lines[i, 9], lines[i, 10], lines[i, 11], lines[i, 12], lines[i, 13],
        lines[i, 14]
      )
    } else {
      s = sprintf(
        "& & & %d & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\ \n",
        lines[i, 4], lines[i, 5], lines[i, 6], lines[i, 7], lines[i, 8], lines[i, 9],
        lines[i, 10], lines[i, 11], lines[i, 12], lines[i, 13], lines[i, 14]
      )
    }
    lines_latex = paste(lines_latex, s)
  }
  return(lines_latex)
}

#===============================================================================

# funcao principal para construcao da tabela ===================================

build_mcse_table <- function(results,
                             indices,
                             n_values = c(10, 25, 50, 100)) {
  
  df <- results %>%
    lapply(function(sc) sc$metrics) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  # para cada indice, filtra, ordena e converte para LaTeX
  latex_lines <- lapply(indices, function(idx) {
    
    df %>%
      filter(n %in% n_values, index == idx) %>%
      select(
        index,
        alpha2,
        true_value,
        n,
        estimator_mean,
        estimator_bias,
        estimator_mse,
        estimator_bias_mcse,
        estimator_mse_mcse,
        corrected_mean,
        corrected_bias,
        corrected_mse,
        corrected_bias_mcse,
        corrected_mse_mcse
      ) %>%
      arrange(alpha2, n) %>%
      converte_cenario_em_latex() %>%
      (function(s) paste(s, "\n \\hline \n "))
  }) %>% do.call(paste, .)
  
  return(latex_lines)
}

#===============================================================================
