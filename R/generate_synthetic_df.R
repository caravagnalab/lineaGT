generate_synthetic_df = function(N_values,
                                 T_values,
                                 K_values,
                                 n_datasets=30,
                                 path=".",
                                 run=T) {

  torch = reticulate::import("torch")
  py_pkg = reticulate::import("pylineaGT")

  files_list = list.files(path=path)

  seeds = sample(0:100, size=n_datasets)

  for (n_df in 1:n_datasets) {
    for (nn in N_values) {
      for (tt in T_values) {
        for (kk in K_values) {
          sim = py_pkg$Simulate(N=as.integer(nn),
                                `T`=as.integer(tt),
                                K=as.integer(kk),
                                seed=as.integer(seeds[n_df]),
                                label=as.integer(n_df))

          if (paste0(sim$sim_id, ".data.Rds") %in% files_list)
            x = readRDS(paste0(sim$sim_id, ".data.Rds"))
          else {
            sim$generate_dataset()
            x = get_simulation_object(sim)
            saveRDS(x, paste0(path, "/", sim$sim_id, ".data.Rds"))
          }

          if (!run)
            next

          labels_true = x$params$z
          cov.df = x$dataset %>%
            tibble::as_tibble() %>%
            dplyr::mutate(labels_true=labels_true,
                          lineage="l1",
                          IS=paste("IS", rownames(.), sep=".")) %>%
            reshape2::melt(id=c("labels_true","IS","lineage"), value.name="coverage", variable.name="timepoints") %>%
            dplyr::mutate(timepoints=str_replace_all(timepoints, "V", "t"))
          k_interval = get_sim_k_interval(x, cov.df)

          x_fit = fit(cov.df=cov.df, k_interval=k_interval, infer_growth=F, infer_phylogenies=F)
          saveRDS(x_fit, paste0(path, "/", sim$sim_id, ".fit.Rds"))

        }
      }
    }
  }
}


get_sim_k_interval = function(x, cov.df) {
  realK = x$settings$K
  k_lower = max(2, min(cov.df %>% long_to_wide_cov() %>% unique() %>% nrow(), realK-5))
  k_upper = min(realK+5, cov.df %>% long_to_wide_cov() %>% unique() %>% nrow())

  return(c(k_lower, k_upper))
}


get_simulation_object = function(sim) {
  x = list()
  x$settings = get_sim_settings(sim$settings)
  x$cov_type = sim$cov_type
  x$dataset = sim$dataset$detach()$numpy()
  x$sim_id = sim$sim_id
  x$params = get_sim_params(sim$params)
  return(x)
}


get_sim_settings = function(settings) {
  sett = list()
  for (name in names(settings))
    tryCatch(
      expr = { sett[[name]] = settings[[name]]$numpy() },
      error = function(e) { sett[[name]] <<- settings[[name]] }
    )

  return(sett)
}


get_sim_params = function(params) {
  pars = list()
  for (name in names(params))
    tryCatch(
      expr = { pars[[name]] = params[[name]]$numpy() },
      error = function(e) { pars[[name]] <<- params[[name]] }
      )

  return(pars)
}


