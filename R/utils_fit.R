# Function to perform a single run of the model
single_fit = function(k, df, run=NULL, steps=500, covariance="diag", lr=0.001, random_state=25) {

  print(paste("RUN", run, "- K =", k))

  x = initialize_object(k, df)
  x = run_inference(x,
                    steps=as.integer(steps),
                    covariance=covariance,
                    lr=as.numeric(lr),
                    random_state=random_state)
  x = classifier(x)

  x$IC = compute_IC(x$py_model)
  x$losses = load_losses(x$py_model)
  x$gradients = load_params_gradients(x$py_model)
  x$n_iter = x$py_model$losses_grad_train$losses %>% length

  return(x)
}


# Function to initialize a python model
# takes as input the long dataframe
initialize_object = function(K, dataset) {

  py_pkg = reticulate::import("pylineaGT")

  lineages = dataset$lineage %>% unique()
  if (!is.null(dataset$timepoints %>% levels())) timepoints = dataset$timepoints %>% levels()
  else timepoints = dataset$timepoints %>% unique()

  df = long_to_wide_input(dataset)
  columns = df %>% dplyr::select(dplyr::starts_with("cov")) %>% colnames()
  IS = df$IS

  py_model = py_pkg$mvnmm$MVNMixtureModel(K=as.integer(K),
                                          data=df %>% dplyr::select(all_of(columns)),
                                          lineages=lineages,
                                          IS=IS,
                                          columns=columns)

  return(get_object(py_model, timepoints=timepoints, lineages=lineages))
}


run_inference = function(x, steps=500, covariance="diag", lr=0.005, random_state=25) {
  x$py_model$fit(steps=as.integer(steps),
                 cov_type=covariance,
                 lr=as.numeric(lr),
                 random_state=as.integer(random_state),
                 convergence=TRUE)
  return(update_params(x))
}


classifier = function(x) {
  x$py_model$classifier()
  return(get_object(x$py_model, timepoints=x$timepoints, lineages=x$lineages))
}

