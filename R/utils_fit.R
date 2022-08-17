# Function to perform a single run of the model
fit_singleK = function(k,
                       cov.df,
                       steps,
                       covariance,
                       hyperparams,
                       lr,
                       p,
                       check_conv,
                       store_params,
                       default_lm,
                       seed_optim,
                       seed,
                       init_seed,
                       timepoints_to_int,
                       py_pkg=NULL) {

  x = initialize_object(K=k,
                        cov.df=cov.df,
                        timepoints_to_int=timepoints_to_int,
                        py_pkg=py_pkg)

  x = run_inference(x,
                    steps=as.integer(steps),
                    lr=as.numeric(lr),

                    check_conv=check_conv,
                    p=as.numeric(p),

                    covariance=covariance,
                    default_lm=default_lm,
                    hyperparams=hyperparams,

                    store_params=store_params,

                    seed_optim=seed_optim,
                    seed=seed,
                    init_seed=init_seed)

  x = classifier(x, timepoints_to_int=timepoints_to_int)

  x$IC = compute_IC(x$py_model)
  x$losses = load_losses(x$py_model)
  x$gradients = load_params_gradients(x$py_model)
  x$n_iter = x$py_model$losses_grad_train$losses %>% length

  x$params$covariance = x$py_model$cov_type

  x$py_model = NULL

  return(x)
}


# Function to initialize a python model
# takes as input the long dataframe
initialize_object = function(K,
                             cov.df,
                             timepoints_to_int=list(),
                             py_pkg=NULL,
                             return_model=FALSE) {
  if (is.null(py_pkg))
    py_pkg = reticulate::import("pylineaGT")

  lineages = cov.df %>% check_lineages()
  timepoints = cov.df %>% check_timepoints()

  df = long_to_wide_cov(cov.df)

  columns = df %>%
    dplyr::select(dplyr::starts_with("cov")) %>%
    colnames()

  IS = df$IS

  py_model = py_pkg$mvnmm$MVNMixtureModel(K=as.integer(K),
                                          data=df %>% dplyr::select(all_of(columns)),
                                          lineages=lineages,
                                          IS=IS,
                                          columns=columns)

  if (return_model) return(py_model)

  return(
    get_object(
      py_model,
      timepoints=timepoints,
      lineages=lineages,
      timepoints_to_int=timepoints_to_int)
    )
}


check_lineages = function(cov.df) {
  if ("lineage" %in% (cov.df %>% colnames()))
    lineages = cov.df$lineage %>% unique()
  else
    lineages = "l.1"

  return(lineages)
}


check_timepoints = function(cov.df) {
  if ("timepoints" %in% (cov.df %>% colnames())) {
    if (!is.null(cov.df$timepoints %>% levels()))
      timepoints = cov.df$timepoints %>% levels()
    else
      timepoints = cov.df$timepoints %>% unique()
  } else
    timepoints = "t.1"

  return(timepoints)
}


run_inference = function(x,
                         steps,
                         lr,

                         covariance,
                         default_lm,
                         hyperparams,

                         check_conv,
                         p,

                         store_params,

                         seed_optim,
                         seed,
                         init_seed) {

  # modify the hyperparams as given in input
  for (hyperpar in names(hyperparams))
    x$py_model$set_hyperparameters(hyperpar, as.numeric(hyperparams[[hyperpar]]))

  x$py_model$fit(steps=as.integer(steps),
                 cov_type=covariance,
                 lr=as.numeric(lr),

                 check_conv=check_conv,
                 p=as.numeric(p),

                 default_lm=default_lm,
                 store_params=store_params,

                 seed_optim=seed_optim,
                 seed=as.integer(seed),
                 init_seed=as.integer(init_seed))

  return(update_params(x))
}


classifier = function(x, timepoints_to_int=list()) {
  x$py_model$classifier()
  return(
    get_object(
      x$py_model,
      timepoints=x$timepoints,
      lineages=x$lineages,
      timepoints_to_int=timepoints_to_int
    )
  )
}

