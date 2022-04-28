# Functions to access elements of a mvnmm object.

get_unique_labels = function(obj) {
  py_model = get_model(obj)
  tryCatch(
    expr = {
      clusters_sort = paste("C_", py_model$params$clusters$detach()$numpy() %>% unique() %>% sort(), sep="")
      return(clusters_sort) },
    error = function(e) return(0:(py_model$params$K - 1)) )
}

get_mean = function(obj) {
  py_model = get_model(obj)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      mean = py_model$params$mean$detach()$numpy()
      colnames(mean) = py_model$dimensions
      rownames(mean) = clusters_sort
      return(mean) },
    error = function(e) return(list()) )
}

get_weights = function(obj) {
  py_model = get_model(obj)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      weights = py_model$params$weights$detach()$numpy()
      names(weights) = clusters_sort
      return(weights) },
    error = function(e) return(list()) )
}

get_sigma = function(obj) {
  py_model = get_model(obj)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      sigma = py_model$params$sigma_vector$detach()$numpy()
      colnames(sigma) = py_model$dimensions
      rownames(sigma) = clusters_sort
      return(sigma) },
    error = function(e) return(list()) )
}

get_covariance_Sigma = function(obj) {
  py_model = get_model(obj)
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
}

get_z_probs = function(obj) {
  py_model = get_model(obj)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      probs = py_model$params$z_probs$detach()$numpy()
      colnames(probs) = clusters_sort
      rownames(probs) = py_model$IS
      return(probs) },
    error = function(e) return(list()) )
}

get_z_assignments = function(obj) {
  py_model = get_model(obj)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      assignments = py_model$params$z_assignments$detach()$numpy()
      colnames(assignments) = clusters_sort
      rownames(assignments) = py_model$IS
      return(assignments) },
    error = function(e) return(list()) )
}

get_labels = function(obj, initial_lab=F) {
  py_model = get_model(obj)
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
}

get_model = function(obj) {
  if ("mvnmm" %in% class(obj)) return(obj$py_model)
  if ("pylineaGT.mvnmm.MVNMixtureModel" %in% class(obj)) return(obj)
}
