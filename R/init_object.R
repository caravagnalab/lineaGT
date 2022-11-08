get_object = function(py_model, timepoints=list(), lineages=list(), timepoints_to_int=list()) {
  x = list()
  x$cov.dataframe = get_python_dataframe(py_model)
  x$params = get_python_params(py_model)

  x$K = py_model$params$K
  x$data.shape = c(py_model$params$N, py_model$params$`T`)

  x$dimensions = py_model$dimensions
  x$timepoints = timepoints
  x$tp.to.int = map_timepoints_int(x, timepoints_to_int=timepoints_to_int)
  x$lineages = lineages

  x$py_model = py_model

  x$color.palette = get_colors(x=x)

  class(x) = "mvnmm"

  return(x)
}


get_python_dataframe = function(py_model) {
  dataset = py_model$dataset$detach()$numpy() %>%
    data.frame() %>%
    dplyr::mutate_all(as.integer)

  try(expr = {
    colnames(dataset) = py_model$dimensions
    dataset$IS = py_model$IS %>% as.character()
    try(expr = {
      labels = suppressMessages( get_labels(py_model) )
      dataset$labels = labels }, silent = T)
  }, silent = T)

  return(dataset %>% wide_to_long_cov())
}


get_python_params = function(py_model) {
  params = list()
  params$mean = suppressMessages( get_mean(py_model) )
  params$weights = suppressMessages( get_weights(py_model) )
  params$sigma = suppressMessages( get_sigma(py_model) )
  params$Sigma = suppressMessages( get_covariance_Sigma(py_model) )
  params$Chol = suppressMessages( get_covariance_Cholesky(py_model) )
  params$probabilites = suppressMessages( get_z_probs(py_model) )
  params$labels = suppressMessages( get_labels(py_model) )
  params$hyperparameters = suppressMessages( get_hyperpar(py_model) )
  params$labels_init = suppressMessages( get_labels(py_model, init=T) )
  return(params)
}


get_hyperpar = function(py_model) {
  hp = list()

  for (hh in names(py_model$hyperparameters))
    hp[[hh]] = py_model$hyperparameters[[hh]]$numpy() %>% as.numeric

  return(hp %>%
           tibble::as_tibble() %>%
           t %>% as.data.frame() %>%
           tibble::rownames_to_column(var="hyperparameter") %>%
           dplyr::rename(value="V1")
         )
}


update_params = function(x) {
  x$params = get_params(py_model=x$py_model)
  x$K = x$py_model$params$K
  x$data.shape = c(x$py_model$params$N, x$py_model$params$`T`)
  return(x)
}



