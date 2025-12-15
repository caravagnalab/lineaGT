generate_synthetic_df = function(N_values,
                                 T_values,
                                 K_values,
                                 path,
                                 n_datasets=30,
                                 likelihood="MVN",
                                 mean_loc=500,
                                 var_loc=110,
                                 var_scale=195,
                                 g_alpha=2,
                                 g_beta=0.05,
                                 steps=1000,
                                 alpha=0.2,
                                 c=3.,
                                 lr=0.005,
                                 filename="",
                                 check_present=TRUE,
                                 default_lm=FALSE,
                                 run=TRUE,
                                 run_NB=FALSE,
                                 py_pkg=NULL) {

  if (!endsWith(path, "/"))
    path = paste0(path, "/")

  if (!dir.exists(path))
    dir.create(path)

  torch = reticulate::import("torch")
  if (is.null(py_pkg)) py_pkg = reticulate::import("pylineaGT")

  alpha_inp = alpha

  seeds = sample(0:min(50,n_datasets), size=n_datasets, replace=FALSE)

  for (n_df in 1:n_datasets) {
    for (nn in N_values) {

      if (!dir.exists(paths=paste0(path, "N", nn)))
        dir.create(paste0(path, "N", nn))

      for (tt in T_values) {
        for (kk in K_values) {

          if (likelihood == "MVN") cli::cli_text("N={nn}, K={kk}, T={tt}, seed={n_df}
                                                 alpha={alpha}, mean_loc={mean_loc},
                                                 var_loc={var_loc}, var_scale={var_scale}")

          else if (likelihood == "NB") cli::cli_text("N={nn}, K={kk}, T={tt}, seed={n_df}
                                                     c={c}, mean_loc={mean_loc},
                                                     g_alpha={g_alpha}, g_beta={g_beta}")

          tmp_name = paste0("N", nn, ".T", tt, ".K", kk)
          subpath = paste0(path, "N", nn, "/", tmp_name, "/")

          if (!dir.exists(paths=subpath))
            dir.create(subpath)

          files_list = list.files(path=subpath)

          if ( check_present &
               (paste0(tmp_name, ".", n_df, ".data.Rds") %in% files_list)  &
               (paste0(tmp_name, ".", n_df, ".fit.Rds") %in% files_list) )
            next

          if ( check_present &
               (paste0(tmp_name, ".", n_df, ".data.Rds") %in% files_list) ) {
            filename = paste0(tmp_name, ".", n_df)
            # filename = paste0(subpath, tmp_name, ".", n_df, ".data.Rds")
            cat(paste0(filename, "\n"))
            x = readRDS(paste0(subpath, filename, ".data.Rds"))
          } else {
            sim = py_pkg$Simulate(N=as.integer(nn),
                                  `T`=as.integer(tt),
                                  K=as.integer(kk),
                                  seed=as.integer(seeds[n_df]),
                                  label=as.character(n_df),
                                  mean_loc=as.integer(mean_loc),
                                  var_loc=as.integer(var_loc),
                                  var_scale=as.integer(var_scale),
                                  g_alpha=as.numeric(g_alpha),
                                  g_beta=as.numeric(g_beta),
                                  alpha=as.numeric(alpha),
                                  c=as.numeric(c),
                                  likelihood=likelihood)

            sim$generate_dataset()
            filename = sim$sim_id

            x = get_simulation_object(sim)

            cat(paste0(subpath, filename, ".data.Rds\n"))
            saveRDS(x, paste0(subpath, filename, ".data.Rds"))
          }

          if (!run & !run_NB) next

          cov.df = x$dataset %>%
            filter_dataset(min_cov=5, min_frac=0)
          k_interval = get_sim_k_interval(x, cov.df)

          if (run) {
            start_time = Sys.time()

            x_fit = fit(cov.df=cov.df,
                        k_interval=k_interval,
                        infer_growth=F,
                        infer_phylogenies=F,
                        covariance="full",
                        check_conv=TRUE,
                        default_lm=default_lm,

                        steps=steps,
                        lr=lr,

                        seed_optim=TRUE,
                        # init_seed=5,
                        sample_id=x$sim_id)

            end_time = Sys.time()

            x_fit$cov.dataframe = tibble::as_tibble(x$dataset) %>%
              dplyr::mutate(coverage=as.integer(coverage)) %>%
              dplyr::inner_join(x_fit$cov.dataframe, by=c("IS","timepoints","lineage","coverage"))

            x_fit$time = end_time - start_time

            print(aricode::NMI(x_fit$cov.dataframe$labels, x_fit$cov.dataframe$labels_true))

            saveRDS(x_fit, paste0(subpath, filename, ".fit.Rds"))
          }

          if (run_NB) {
            start_time = Sys.time()

            cov.df_wide = cov.df %>% tidyr::unite("covariate", c("lineage","timepoints"), sep="_") %>%
              tidyr::pivot_wider(id_cols=c("labels_true", "IS"),
                                 names_from="covariate", values_from="coverage")

            # input_df = cov.df_wide %>% select(-labels_true, -IS)
            # Ks = replicate(10, {
            #   NB.MClust::NB.MClust(Count=input_df, K=k_interval, iteration=500)$K
            # })
            #
            # best_k = which.max(table(Ks)) %>% names() %>% as.integer()
            # x_fit.NB = NB.MClust::NB.MClust(Count=input_df, K=best_k, iteration=10000)
            #
            # cov.df_wide$labels = x_fit.NB$cluster
            # x_fit.NB$cov.dataframe = cov.df_wide %>%
            #   tidyr::pivot_longer(cols=-c("labels", "labels_true", "IS"), values_to="coverage") %>%
            #   tidyr::separate(name, into=c("lineage", "timepoints"))
            #
            # end_time = Sys.time()
            #
            # x_fit.NB$time = end_time - start_time
            #
            # print(aricode::NMI(x_fit.NB$cov.dataframe$labels, x_fit.NB$cov.dataframe$labels_true))
            # saveRDS(x_fit, paste0(subpath, filename, ".fit_NB.Rds"))
          }

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


