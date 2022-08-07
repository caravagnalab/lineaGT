# Function to perform a single run of the model
fit_singleK = function(k,
                       cov.df,
                       steps=500,
                       covariance="full",
                       hyperparameters=list(),
                       lr=0.001,
                       p=1,
                       convergence=TRUE,
                       store_params=FALSE,
                       initializ=TRUE,
                       default_constr=TRUE,
                       sigma_constr_pars=list("slope"=0.09804862, "intercept"=22.09327233),
                       seed=5,
                       timepoints_to_int=list(),
                       py_pkg=NULL) {

  x = initialize_object(K=k,
                        default_constr=default_constr,
                        sigma_constr_pars=sigma_constr_pars,
                        cov.df=cov.df,
                        py_pkg=py_pkg,
                        timepoints_to_int=timepoints_to_int)

  x = run_inference(x,
                    steps=as.integer(steps),
                    covariance=covariance,
                    lr=as.numeric(lr),
                    p=as.numeric(p),
                    hyperparameters=hyperparameters,
                    convergence=convergence,
                    store_params=store_params,
                    initializ=initializ,
                    seed=seed)

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
                             default_constr=TRUE,
                             sigma_constr_pars=list("slope"=0.09804862, "intercept"=22.09327233),
                             timepoints_to_int=list(),
                             py_pkg=NULL) {
  if (is.null(py_pkg))
    py_pkg = reticulate::import("pylineaGT")

  sigma_constr_pars = reticulate::py_dict(keys=names(sigma_constr_pars),
                                          values=as.numeric(sigma_constr_pars))

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
                                          columns=columns,
                                          default_init=default_constr)



  if (default_constr)
    py_model$set_sigma_constraints(slope=sigma_constr_pars["slope"], intercept=sigma_constr_pars["intercept"])

  return(get_object(py_model, timepoints=timepoints, lineages=lineages, timepoints_to_int=timepoints_to_int))
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
                         steps=500,
                         covariance="full",
                         hyperparameters=list(),
                         lr=0.005,
                         p=1,
                         convergence=TRUE,
                         initializ=TRUE,
                         store_params=FALSE,
                         seed=5) {

  # modify the hyperparameters as given in input
  for (hyperpar in names(hyperparameters))
    x$py_model$set_hyperparameters(hyperpar, as.numeric(hyperparameters[[hyperpar]]))

  x$py_model$fit(steps=as.integer(steps),
                 cov_type=covariance,
                 lr=as.numeric(lr),
                 p=as.numeric(p),
                 convergence=convergence,
                 store_params=store_params,
                 initializ=initializ,
                 seed=as.integer(seed))

  return(update_params(x))
}


classifier = function(x, timepoints_to_int=list()) {
  x$py_model$classifier()
  return(get_object(x$py_model, timepoints=x$timepoints, lineages=x$lineages, timepoints_to_int=timepoints_to_int))
}

