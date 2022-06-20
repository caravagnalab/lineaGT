# Function to perform a single run of the model
fit_singleK = function(k,
                       cov.df,
                       steps=500,
                       covariance="diag",
                       lr=0.001,
                       p=0.01,
                       convergence=TRUE,
                       random_state=25,
                       py_pkg=NULL) {

  x = initialize_object(k, cov.df, py_pkg)
  x = run_inference(x,
                    steps=as.integer(steps),
                    covariance=covariance,
                    lr=as.numeric(lr),
                    p=as.numeric(p),
                    convergence=convergence,
                    random_state=random_state)
  x = classifier(x)

  x$IC = compute_IC(x$py_model)
  x$losses = load_losses(x$py_model)
  x$gradients = load_params_gradients(x$py_model)
  x$n_iter = x$py_model$losses_grad_train$losses %>% length

  x$py_model = NULL

  return(x)
}


# Function to initialize a python model
# takes as input the long dataframe
initialize_object = function(K,
                             cov.df,
                             py_pkg=NULL) {
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

  return(get_object(py_model, timepoints=timepoints, lineages=lineages))
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


run_inference = function(x, steps=500, covariance="diag", lr=0.005,
                         p=0.01, convergence=TRUE, random_state=25) {
  x$py_model$fit(steps=as.integer(steps),
                 cov_type=covariance,
                 lr=as.numeric(lr),
                 p=as.numeric(p),
                 convergence=convergence,
                 random_state=as.integer(random_state))
  return(update_params(x))
}


classifier = function(x) {
  x$py_model$classifier()
  return(get_object(x$py_model, timepoints=x$timepoints, lineages=x$lineages))
}

