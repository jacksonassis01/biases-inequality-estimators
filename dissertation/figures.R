#===============================================================================
# figures.R
# Geracao de figuras de Vies e EQM a partir de "results"
#===============================================================================

library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

# funcao auxiliar de-para nos labels dos indices ===============================
pretty_index_name <- function(idx) {
  switch(
    idx,
    "TheilT"      = "Theil-T",
    "TheilL"      = "Theil-L",
    "Atkinson1"   = "Atkinson (ε = 1)",
    "AtkinsonInf" = "Atkinson (ε → ∞)",
    "VMR"         = "VMR",
    idx
  )
}

#===============================================================================

# funcoes internas de plotagem =================================================
build_row_plot <- function(data_bias, data_mse, indexselected, xvar) {
  g_bias <- data_bias %>%
    ggplot() +
    aes_string(x = xvar, y = "bias", linetype = "kind") +
    geom_line(color = "gray32") +
    geom_point(color = "gray32") +
    scale_linetype_manual(
      NULL,
      values = c(
        "Estimador" = "dashed",
        "Estimador Corrigido" = "solid"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      plot.margin = margin(
        t = 0,
        b = 0,
        l = 0,
        r = 5
      ),
      panel.grid.major = element_line(color = "gray96"),
      panel.grid.minor = element_line(color = "gray96")
    )
  
  g_mse <- data_mse %>%
    ggplot() +
    aes_string(x = xvar, y = "mse", linetype = "kind") +
    geom_line(color = "gray32") +
    geom_point(color = "gray32") +
    scale_linetype_manual(
      NULL,
      values = c(
        "Estimador" = "dashed",
        "Estimador Corrigido" = "solid"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      plot.margin = margin(
        t = 0,
        b = 0,
        l = 0,
        r = 5
      ),
      panel.grid.major = element_line(color = "gray96"),
      panel.grid.minor = element_line(color = "gray96")
    )
  
  g_title <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = pretty_index_name(indexselected),
      angle = 90,
      size = 3.5,
      hjust = 0.5,
      vjust = 0.5,
      color = "gray20",
      fontface = "bold"
    ) +
    theme_void() +
    theme(plot.margin = margin(r = 2, l = 2))
  
  g_title + g_bias + g_mse +
    plot_layout(widths = c(0.05, 1, 1), design = "ABC")
}

assemble_final_figure <- function(rows) {
  header_bias <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "Viés",
      size = 4,
      fontface = "bold",
      color = "gray20"
    ) +
    theme_void()
  header_eqm <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "EQM",
      size = 4,
      fontface = "bold",
      color = "gray20"
    ) +
    theme_void()
  header_blank <- ggplot() + theme_void()
  
  header <- header_blank + header_bias + header_eqm +
    plot_layout(widths = c(0.05, 1, 1))
  
  rows_combined <- wrap_plots(plotlist = rows, ncol = 1)
  
  wrap_plots(
    header,
    rows_combined,
    ncol = 1,
    heights = c(0.05, 1),
    guides = "collect"
  ) &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      axis.title = element_blank()
    )
}

#===============================================================================

# conversao de results em dataframe ============================================
results_to_df <- function(results) {
  do.call(rbind, lapply(results, function(obj)
    obj$metrics)) %>%
    as.data.frame()
}

#===============================================================================

# funcoes principais de plotagem ===============================================

# plot por n, fixando alpha1 e alpha2
plot_indices_by_n_from_results <- function(results,
                                           alpha1_selected,
                                           alpha2_selected,
                                           indices) {
  df_all <- results_to_df(results)
  n_range <- unique(df_all$n)
  
  rows <- lapply(indices, function(idx) {
    df_bias <- df_all %>%
      filter(index == idx,
             alpha1 == alpha1_selected,
             alpha2 == alpha2_selected,
             n %in% n_range) %>%
      select(index, n, alpha1, alpha2, estimator_bias, corrected_bias) %>%
      melt(
        id.vars = c("index", "n", "alpha1", "alpha2"),
        variable.name = "type",
        value.name = "bias"
      ) %>%
      mutate(kind = ifelse(
        type == "estimator_bias",
        "Estimador",
        "Estimador Corrigido"
      ))
    
    df_mse <- df_all %>%
      filter(index == idx,
             alpha1 == alpha1_selected,
             alpha2 == alpha2_selected,
             n %in% n_range) %>%
      select(index, n, alpha1, alpha2, estimator_mse, corrected_mse) %>%
      melt(
        id.vars = c("index", "n", "alpha1", "alpha2"),
        variable.name = "type",
        value.name = "mse"
      ) %>%
      mutate(kind = ifelse(type == "estimator_mse", "Estimador", "Estimador Corrigido"))
    
    build_row_plot(df_bias, df_mse, idx, "n")
  })
  
  assemble_final_figure(rows)
}

# plot por alpha2, fixando n e alpha1
plot_indices_by_alpha2_from_results <- function(results,
                                                n_selected,
                                                alpha1_selected,
                                                indices) {
  df_all <- results_to_df(results)
  
  rows <- lapply(indices, function(idx) {
    df_bias <- df_all %>%
      filter(index == idx, n == n_selected, alpha1 == alpha1_selected) %>%
      select(index, n, alpha1, alpha2, estimator_bias, corrected_bias) %>%
      melt(
        id.vars = c("index", "n", "alpha1", "alpha2"),
        variable.name = "type",
        value.name = "bias"
      ) %>%
      mutate(kind = ifelse(
        type == "estimator_bias",
        "Estimador",
        "Estimador Corrigido"
      ))
    
    df_mse <- df_all %>%
      filter(index == idx, n == n_selected, alpha1 == alpha1_selected) %>%
      select(index, n, alpha1, alpha2, estimator_mse, corrected_mse) %>%
      melt(
        id.vars = c("index", "n", "alpha1", "alpha2"),
        variable.name = "type",
        value.name = "mse"
      ) %>%
      mutate(kind = ifelse(type == "estimator_mse", "Estimador", "Estimador Corrigido"))
    
    build_row_plot(df_bias, df_mse, idx, "alpha2")
    
  })
  
  assemble_final_figure(rows)
}

#===============================================================================
