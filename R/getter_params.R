get_params = function(x=NULL, py_model=NULL) {
  if (!is.null(x)) return(x$params)

  if(!is.null(py_model)) return(get_python_params(py_model))
}



get_mean = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      mean = py_model$params$mean$detach()$numpy()
      colnames(mean) = py_model$dimensions
      rownames(mean) = clusters_sort
      print("Done")
      return(mean) },
    error = function(e) return(list()) )

  try(expr = { mean = x$params$mean; if (!purrr::is_empty(mean)) return(mean) }, silent = T)
}



get_weights = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      weights = py_model$params$weights$detach()$numpy()
      names(weights) = clusters_sort
      return(weights) },
    error = function(e) return(list()) )

  try(expr = { weights = x$params$weights; if (!purrr::is_empty(weights)) return(weights) }, silent = T)
}



get_sigma = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      sigma = py_model$params$sigma_vector$detach()$numpy()
      colnames(sigma) = py_model$dimensions
      rownames(sigma) = clusters_sort
      return(sigma) },
    error = function(e) return(list()) )

  try(expr = { sigma = x$params$sigma; if (!purrr::is_empty(sigma)) return(sigma) }, silent = T)
}



get_covariance_Sigma = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      covar = list()
      for (k in 0:(py_model$params$K-1)) {
        name = paste("C_", k, sep="")
        covar[[name]] = py_model$params$sigma[k]$detach()$numpy()
        covar[[name]] = covar[[name]] %*% t(covar[[name]])
        colnames(covar[[name]]) = rownames(covar[[name]]) = py_model$dimensions
      }
      return(covar) },
    error = function(e) return(list()) )

  try(expr = { Sigma = x$params$Sigma; if (!purrr::is_empty(Sigma)) return(Sigma) }, silent = T)
}



get_z_probs = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      probs = py_model$params$z_probs$detach()$numpy()
      colnames(probs) = clusters_sort
      rownames(probs) = py_model$IS
      return(probs) },
    error = function(e) return(list()) )

  try(expr = { z_probs = x$params$probabilites; if (!purrr::is_empty(z_probs)) return(z_probs) }, silent = T)
}



get_z_assignments = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      assignments = py_model$params$z_assignments$detach()$numpy()
      colnames(assignments) = clusters_sort
      rownames(assignments) = py_model$IS
      return(assignments) },
    error = function(e) return(list()) )

  try(expr = { assignm = x$params$assignments; if (!purrr::is_empty(assignm)) return(assignm) }, silent = T)
}



get_labels = function(x, initial_lab=F) {
  py_model = get_model(x)
  if (!initial_lab) {
    tryCatch(
      expr = {
        clusters_sort = get_unique_labels(py_model)
        labels = factor(paste("C_", py_model$params$clusters$detach()$numpy(), sep=""), levels=clusters_sort)
        return(labels) },
      error = function(e) return(list()) ) }

  else {
    tryCatch(
      expr = {
        clusters_sort = get_unique_labels(py_model)
        labels = factor(paste("C_", py_model$init_params$clusters$detach()$numpy(), sep=""), levels=clusters_sort)
        return(labels) },
      error = function(e) return(list()) ) }

  try(expr = { labs = x$params$labels; if (!purrr::is_empty(labs)) return(labs) }, silent = T)
}



get_unique_labels = function(x) {
  try(expr = { labels = x$params$labels %>% levels(); if (!purrr::is_empty(labels)) return(labels) }, silent = T)

  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = paste("C_", py_model$params$clusters$detach()$numpy() %>% unique() %>% sort(), sep="")
      return(clusters_sort) },
    error = function(e) return(0:(py_model$params$K - 1)) )
}



get_unique_viber_labels = function(x) {
  labels = x$vaf.dataframe$labels_mut %>% unique()
  return(labels)
}
