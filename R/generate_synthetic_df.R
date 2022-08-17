generate_synthetic_df = function(N_values,
                                 T_values,
                                 K_values,
                                 n_datasets=30,
                                 var_loc=118,
                                 var_scale=130,
                                 mean_loc=500,
                                 mean_scale=1000,
                                 path=".",
                                 filename="",
                                 run=T) {

  torch = reticulate::import("torch")
  py_pkg = reticulate::import("pylineaGT")

  files_list = list.files(path=path)

  seeds = sample(0:min(50,n_datasets), size=n_datasets, replace=FALSE)

  for (n_df in 1:n_datasets) {
    for (nn in N_values) {
      for (tt in T_values) {
        for (kk in K_values) {
          sim = py_pkg$Simulate(N=as.integer(nn),
                                `T`=as.integer(tt),
                                K=as.integer(kk),
                                seed=as.integer(seeds[n_df]),
                                label=as.integer(n_df),
                                var_loc=as.integer(var_loc),
                                var_scale=as.integer(var_scale),
                                mean_loc=as.integer(mean_loc),
                                mean_scale=as.integer(mean_scale))


          if (filename == "") filename = sim$sim_id

          if (paste0(filename, ".data.Rds") %in% files_list) x = readRDS(paste0(path, filename, ".data.Rds"))
          else {
            sim$generate_dataset()
            x = get_simulation_object(sim)

            print(paste0(path, filename, ".data.Rds"))
            saveRDS(x, paste0(path, filename, ".data.Rds"))
          }

          if (!run) next

          cov.df = x$dataset %>%
            filter_dataset(min_cov=5, min_frac=0)

          k_interval = get_sim_k_interval(x, cov.df)

          x_fit = fit(cov.df=cov.df,
                      k_interval=k_interval,
                      infer_growth=F,
                      infer_phylogenies=F,
                      default_lm=TRUE,
                      seed_optim=TRUE,
                      sample_id=x$sim_id)

          x_fit$cov.dataframe = tibble::as_tibble(x$dataset) %>%
            dplyr::mutate(coverage=as.integer(coverage)) %>%
            dplyr::inner_join(x_fit$cov.dataframe, by=c("IS","timepoints","lineage","coverage"))

          saveRDS(x_fit, paste0(path, filename, ".fit.Rds"))
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
  x$params = get_sim_params(sim$params)
  labels_true = x$params$z
  x$settings = get_sim_settings(sim$settings)
  x$cov_type = sim$cov_type
  x$dataset = get_sim_dataset(sim, labels_true)
  x$sim_id = sim$sim_id

  return(x)
}


get_sim_dataset = function(sim, labels_true){
  dataset = sim$dataset$detach()$numpy()
  df = dataset %>%
    tibble::as_tibble() %>%
    dplyr::mutate(labels_true=labels_true,
                  lineage="l1",
                  IS=paste("IS", rownames(.), sep=".")) %>%
    reshape2::melt(id=c("labels_true","IS","lineage"), value.name="coverage", variable.name="timepoints") %>%
    dplyr::mutate(timepoints=str_replace_all(timepoints, "V", "t"))
  return(df)
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

  pars$z = paste0("C", pars$z)

  return(pars)
}


