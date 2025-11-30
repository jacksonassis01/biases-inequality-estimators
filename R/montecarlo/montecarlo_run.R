# R/montecarlo/montecarlo_run.R

# Funções para executar simulações Monte Carlo para mistura Gama

run_montecarlo <- function(R,
                           n_range,
                           alpha1_values,
                           alpha2_values,
                           pii,
                           beta,
                           indices,
                           mle_fun,
                           rng_fun,
                           save_tmp = FALSE,
                           tmp_dir = "tmp") {
  stopifnot(is.list(indices), length(indices) > 0)
  
  if (save_tmp &&
      !dir.exists(tmp_dir))
    dir.create(tmp_dir, recursive = TRUE)
  
  # Lista de parâmetros verdadeiros (todos os cenários)
  
  params_true_all <- list()
  for (n in n_range) {
    for (a1 in alpha1_values) {
      for (a2 in alpha2_values) {
        nome <- sprintf("alpha1_%.2f_alpha2_%.2f_n_%d", a1, a2, n)
        params_true_all[[nome]] <- list(
          n     = n,
          pii   = pii,
          alpha = c(a1, a2),
          beta  = beta
        )
      }
    }
  }
  
  # Wrapper para gerar amostras
  
  rpop <- function(fun, args_list)
    do.call(fun, args_list)
  
  # Loop sobre cenários
  
  results <- lapply(names(params_true_all), function(nome) {
    params_true <- params_true_all[[nome]]
    n <- params_true$n
    a1 <- params_true$alpha[1]
    a2 <- params_true$alpha[2]
    
    # Valores verdadeiros dos índices
    theta_true_list <- lapply(indices, function(idx) {
      args <- params_true[names(params_true) %in% names(formals(idx$true_func))]
      do.call(idx$true_func, args)
    })
    
    # Gerar amostras
    samples <- replicate(R, rpop(
      rng_fun,
      list(
        n = n,
        pii = params_true$pii,
        alpha = params_true$alpha,
        beta = params_true$beta
      )
    ), simplify = FALSE)
    
    # Rodar estimadores e corrigir bias
    sims <- lapply(samples, function(sample) {
      emout <- mle_fun(sample)
      est_params <- list(
        pii   = c(emout$pi1, emout$pi2),
        alpha = c(emout$alpha1, emout$alpha2),
        beta  = emout$beta,
        n     = n
      )
      
      index_results <- lapply(names(indices), function(name) {
        idx <- indices[[name]]
        theta_hat <- idx$est_func(sample)
        bias_est <- do.call(idx$bias_func, est_params)
        theta_corr <- theta_hat - bias_est
        list(
          index = name,
          theta_hat = theta_hat,
          theta_corr = theta_corr,
          bias_est = bias_est
        )
      })
      
      list(
        z = sample,
        estimated_params = est_params,
        estimated_indices = index_results
      )
    })
    
    result <- list(
      params_true = params_true,
      theta_true  = theta_true_list,
      simulations = sims
    )
    
    # Salvar temporariamente se solicitado
    if (save_tmp) {
      filename <- file.path(tmp_dir, paste0(nome, "_R", R, ".rds"))
      saveRDS(result, filename)
    }
    
    result
  })
  
  names(results) <- names(params_true_all)
  return(results)
}

run_full_experiment <- function(alpha1_values = c(0.5),
                                alpha2_values = c(0.5),
                                n_values = c(10, 50),
                                Nrep = 1000,
                                pii = c(0.6, 0.4),
                                beta = 1,
                                indices = default_indices_list,
                                mle_fun = estimate_gamma_mixture_ls,
                                rng_fun = rMGamma,
                                verbose = TRUE,
                                save_tmp = FALSE,
                                tmp_dir = "tmp") {
  if (save_tmp &&
      !dir.exists(tmp_dir))
    dir.create(tmp_dir, recursive = TRUE)
  
  total_runs <- length(alpha1_values) * length(alpha2_values) * length(n_values)
  run_counter <- 1
  
  all_results <- list()
  
  for (a1 in alpha1_values) {
    for (a2 in alpha2_values) {
      for (n in n_values) {
        start_time <- Sys.time()
        if (verbose)
          cat(
            sprintf(
              "[%d/%d] Running simulation for alpha1 = %.2f, alpha2 = %.2f, n = %d ... ",
              run_counter,
              total_runs,
              a1,
              a2,
              n
            )
          )
        
        sim_out <- run_montecarlo(
          R = Nrep,
          n_range = n,
          alpha1_values = a1,
          alpha2_values = a2,
          pii = pii,
          beta = beta,
          indices = indices,
          mle_fun = mle_fun,
          rng_fun = rng_fun,
          save_tmp = save_tmp,
          tmp_dir = tmp_dir
        )
        
        key <- sprintf("alpha1_%.2f_alpha2_%.2f_n_%d", a1, a2, n)
        all_results[[key]] <- sim_out[[1]]  # sim_out sempre tem apenas 1 elemento por chamada
        
        elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
        if (verbose)
          cat(sprintf("done (%.2f sec)\n", elapsed))
        run_counter <- run_counter + 1
      }
    }
  }
  
  if (verbose)
    cat("All simulations have been completed!\n")
  all_results
}

append_metrics <- function(results) {
  compute_stats <- function(vec, true_val, Nrep) {
    diff <- vec - true_val
    diff2 <- diff ^ 2
    list(
      mean_value = mean(vec),
      bias_value = mean(diff),
      mse_value = mean(diff2),
      bias_mcse = sd(diff) / sqrt(Nrep),
      mse_mcse = sd(diff2) / sqrt(Nrep)
    )
  }
  
  compute_metrics <- function(obj) {
    n <- obj$params_true$n
    alpha1 <- obj$params_true$alpha[1]
    alpha2 <- obj$params_true$alpha[2]
    true_vals <- obj$theta_true
    
    sims_list <- lapply(seq_along(obj$simulations), function(i) {
      sim <- obj$simulations[[i]]
      idx_list <- lapply(sim$estimated_indices, function(idx) {
        data.frame(
          sim_id = i,
          index = idx$index,
          theta_hat = idx$theta_hat,
          theta_corr = idx$theta_corr,
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, idx_list)
    })
    
    sim_df <- do.call(rbind, sims_list)
    Nrep <- length(unique(sim_df$sim_id))
    
    metrics <- do.call(rbind, lapply(split(sim_df, sim_df$index), function(df) {
      idx <- unique(df$index)
      true_value <- true_vals[[idx]]
      stats_estimator <- compute_stats(df$theta_hat, true_value, Nrep)
      stats_corrected <- compute_stats(df$theta_corr, true_value, Nrep)
      
      data.frame(
        index = idx,
        n = n,
        alpha1 = alpha1,
        alpha2 = alpha2,
        estimator_mean = stats_estimator$mean_value,
        estimator_bias = stats_estimator$bias_value,
        estimator_mse = stats_estimator$mse_value,
        estimator_bias_mcse = stats_estimator$bias_mcse,
        estimator_mse_mcse = stats_estimator$mse_mcse,
        corrected_mean = stats_corrected$mean_value,
        corrected_bias = stats_corrected$bias_value,
        corrected_mse = stats_corrected$mse_value,
        corrected_bias_mcse = stats_corrected$bias_mcse,
        corrected_mse_mcse = stats_corrected$mse_mcse,
        true_value = true_value,
        stringsAsFactors = FALSE
      )
    }))
    
    rownames(metrics) <- NULL
    
    obj$metrics <- metrics
    obj
  }
  
  lapply(results, compute_metrics)
}
