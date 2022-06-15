get_params = function(x=NULL, py_model=NULL) {
  if (!is.null(x)) return(x$params)

  if(!is.null(py_model)) return(get_python_params(py_model))
}


#' Extract the estimated mean parameters.
#'
#' @description Returns a dataframe \code{KxT} with the estimated mean paramaters \code{mu_kt} per
#' clone \code{k} and dimension \code{t}.
#'
#' @param x a mvnmm object.
#' @return dataframe of mean parameters.
#'
#' @examples
#' if (FALSE) get_mean(x)
#'
#' @export get_mean

get_mean = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model)
      mean = py_model$params$mean$detach()$numpy()
      colnames(mean) = py_model$dimensions
      rownames(mean) = clusters_sort
      return(mean) },
    error = function(e) return(list()) )

  try(expr = { mean = x$params$mean; if (!purrr::is_empty(mean)) return(mean) }, silent = T)
}


#' Extract the estimated mixing proportions.
#'
#' @description Returns a list of dimension \code{K} with the estimated mixing proportions.
#'
#' @param x a mvnmm object.
#' @return list of estimated mixing proportions.
#'
#' @examples
#' if (FALSE) get_weights(x)
#'
#' @export get_weights

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


#' Extract the estimated variance parameters.
#'
#' @description Returns a dataframe \code{KxT} with the estimated variance paramaters \code{sigma_kt} per
#' clone \code{k} and dimension \code{t}.
#'
#' @param x a mvnmm object.
#' @return dataframe of variance parameters.
#'
#' @examples
#' if (FALSE) get_sigma(x)
#'
#' @export get_sigma

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


#' Extract the estimated covariance matrices.
#'
#' @description Returns a list with \code{K} dataframes, each of dimension \code{TxT}, corresponding
#' to the covariance matrices estimated for each clone \code{k}
#'
#' @param x a mvnmm object.
#' @return list of the estimated covariance matrices.
#'
#' @examples
#' if (FALSE) get_covariance_Sigma(x)
#'
#' @export get_covariance_Sigma

get_covariance_Sigma = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      covar = list()
      for (k in 0:(x$K-1)) {
        name = paste("C_", k, sep="")
        covar[[name]] = py_model$params$sigma[k]$detach()$numpy()
        covar[[name]] = covar[[name]] %*% t(covar[[name]])
        colnames(covar[[name]]) = rownames(covar[[name]]) = py_model$dimensions
      }
      return(covar) },
    error = function(e) return(list()) )

  try(expr = { Sigma = x$params$Sigma; if (!purrr::is_empty(Sigma)) return(Sigma) }, silent = T)
}


#' Extract the estimated posterior probabilities.
#'
#' @description Returns a dataframe of shape \code{NxK} with the posterior distribution \code{p(k|n)}
#' for each observation \code{n} to belong to cluster \code{k}.
#'
#' @param x a mvnmm object.
#' @return dataframe of posterior distributions.
#'
#' @examples
#' if (FALSE) get_z_probs(x)
#'
#' @export get_z_probs

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


#' Extract the observations labels.
#'
#' @description Returns a list with \code{N} elements, corresponding to the labels for each
#' observation.
#'
#' @param x a mvnmm object.
#' @return list of observations labels.
#'
#' @examples
#' if (FALSE) get_labels(x)
#'
#' @export get_labels

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


#' Extract the list of unique observations labels.
#'
#' @description Returns a list with \code{K} elements, corresponding to the unique labels.
#'
#' @param x a mvnmm object.
#' @return list of unique labels.
#'
#' @examples
#' if (FALSE) get_unique_labels(x)
#'
#' @export get_unique_labels

get_unique_labels = function(x) {
  try(expr = { labels = x$params$labels %>% levels(); if (!purrr::is_empty(labels)) return(labels) }, silent = T)

  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = paste("C_", py_model$params$clusters$detach()$numpy() %>% unique() %>% sort(), sep="")
      return(clusters_sort) },
    error = function(e) return(0:(py_model$params$K - 1)) )
}


#' Retrieve the list of unique labels of mutation clusters.
#'
#' @description Function to retrieve the list of unique labels of mutations clusters,
#' of the form \code{C_c1.Cm1}, where \code{c1} is the clone identifier and \code{m1}
#' is the subclone identifier.
#'
#' @param x a mvnmm object.
#' @param clusters a vector-like variable, with the identifiers of the clones we want
#' to retrieve the subclone labels from. If empty, all the labels will be returned.
#' @return vector of mutations labels.
#'
#' @examples
#' if(FALSE) get_unique_muts_labels(x, c("C_0"))
#'
#' @export get_unique_muts_labels

get_unique_muts_labels = function(x, clusters=c(), label="") {
  if (purrr::is_empty(clusters)) return(get_all_unique_muts_labels(x, label=label))
  return(
    x %>%
      get_vaf_dataframe(label=label) %>%
      dplyr::filter(labels %in% clusters) %>%
      dplyr::pull(labels_mut) %>%
      unique()
  )
}


get_all_unique_muts_labels = function(x, label="") {
  labels = x %>%
    get_vaf_dataframe(label=label) %>%
    dplyr::pull(labels_mut) %>%
    unique()
  return(labels)
}

