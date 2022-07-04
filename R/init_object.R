get_object = function(py_model, timepoints=list(), lineages=list()) {
  x = list()
  x$cov.dataframe = get_python_dataframe(py_model)
  x$params = get_python_params(py_model)

  x$K = py_model$params$K
  x$N = py_model$params$N
  x$`T` = py_model$params$`T`

  x$dimensions = py_model$dimensions
  x$timepoints = timepoints
  x$tp_to_int = map_timepoints_int(x)
  x$lineages = lineages

  x$py_model = py_model

  x$color_palette = get_colors(x=x)

  class(x) = "mvnmm"

  return(x)
}


get_python_dataframe = function(py_model) {
  dataset = py_model$dataset$detach()$numpy() %>%
    data.frame() %>%
    dplyr::mutate_all(as.integer)

  try(expr = {
    colnames(dataset) = py_model$dimensions
    dataset$IS = py_model$IS
    try(expr = {
      labels = get_labels(py_model)
      dataset$labels = labels }, silent = T)

    try(expr = {
      labels_init = get_labels(py_model, initial_lab=T)
      dataset$labels_init = labels_init }, silent = T)
  }, silent = T)

  return(dataset %>% wide_to_long_cov())
}


get_python_params = function(py_model, train=FALSE) {
  params = list()
  if (!train) {
    params$mean = get_mean(py_model)
    params$weights = get_weights(py_model)
    params$sigma = get_sigma(py_model)
    params$Sigma = get_covariance_Sigma(py_model)
    params$probabilites = get_z_probs(py_model)
    params$labels = get_labels(py_model)
    params$hyperparameters = py_model$hyperparameters
  } else {

  }
  return(params)
}


update_params = function(x) {
  x$params = get_params(py_model=x$py_model)
  x$params_train = get_params(py_model=x$py_model, train=T)
  x$K = x$py_model$params$K
  x$N = x$py_model$params$N
  x$`T` = x$py_model$params$`T`
  return(x)
}



