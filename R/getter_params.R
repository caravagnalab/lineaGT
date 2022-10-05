get_default_hyperparameters = function() {
  py_pkg = reticulate::import("pylineaGT")
  x = initialize_object(K=1,
                        cov.df=data.frame(coverage=1:10, IS=1:10, timepoints=1, lineage=1),
                        py_pkg=py_pkg)

  hyper = x$params$hyperparameters

  return(hyper %>% dplyr::filter(!hyperparameter %in%c("mean_scale","mean_loc")))
}


# function to get parameters from either the Python model or the mvnmm object
get_params = function(x=NULL, py_model=NULL) {
  if (!is.null(x) && !train) return( x$params )
  if(!is.null(py_model)) return( get_python_params(py_model) )
  cli::format_error("Not able to return any 'params' object!")
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
  tryCatch(  # case 1 -> we want to get the params from the Python object
    expr = {
      clusters_sort = get_unique_labels(py_model)
      mean = py_model$params$mean$detach()$numpy()
      colnames(mean) = py_model$dimensions
      rownames(mean) = clusters_sort
      return(mean) },
    error = function(e) { return(list()) } )

  check_params(x)
  if ("mean" %in% names(x$params)) return(x$params$mean)
  cli::cli_alert_warning("No mean stored in the object.")
  return(list())
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

  check_params(x)
  if ("weights" %in% names(x$params)) return(x$params$weights)

  cli::cli_alert_warning("No weights stored in the object.")
  return(list())
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

  check_params(x)
  if ("sigma" %in% names(x$params)) return(x$params$sigma)

  cli::cli_alert_warning("No variance 'sigma' stored in the object.")
  return(list())
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
        name = paste("C", k, sep="")
        covar[[name]] = py_model$params$sigma[k]$detach()$numpy()
        covar[[name]] = covar[[name]] %*% t(covar[[name]])
        colnames(covar[[name]]) = rownames(covar[[name]]) = py_model$dimensions
      }
      return(covar) },
    error = function(e) return(list()) )

  check_params(x)
  if ("Sigma" %in% names(x$params)) return(x$params$Sigma)

  cli::cli_alert_warning("No covariance 'Sigma' stored in the object.")
  return(list())
}


#' Extract the estimated Cholesky matrices, used to factorise the covariance matrix.
#'
#' @description Returns a list with \code{K} dataframes, each of dimension \code{TxT}, corresponding
#' to the covariance matrices estimated for each clone \code{k}
#'
#' @param x a mvnmm object.
#' @return list of the estimated covariance matrices.
#'
#' @examples
#' if (FALSE) get_covariance_Cholesky(x)
#'
#' @export get_covariance_Cholesky

get_covariance_Cholesky = function(x) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      if (py_model$cov_type=="diag") {
        chol = py_model$params$sigma_chol$detach()$numpy()
        colnames(chol) = rownames(chol) = py_model$dimensions
        return(chol)
      }

      chol = list()
      for (k in 0:(x$K-1)) {
        name = paste("C", k, sep="")
        chol[[name]] = py_model$params$sigma_chol[k]$detach()$numpy()
        colnames(chol[[name]]) = rownames(chol[[name]]) = py_model$dimensions
      }
      return(chol) },
    error = function(e) return(list()) )

  check_params(x)
  if ("Chol" %in% names(x$params)) return(x$params$Chol)

  cli::cli_alert_warning("No Cholesky factorization matrix 'Chol' stored in the object.")
  return(list())
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

  check_params(x)
  if ("probabilites" %in% names(x$params)) return(x$params$probabilites)

  cli::cli_alert_warning("No posterior probabilities stored in the object.")
  return(list())
}


#' Get the number of ISs per cluster.
#'
#' @param x the fitted object
#' @param highlight the clusters to retrieve
#'
#' @return an array with as names the clusters in \code{highlight} and as values the number of ISs assigned
#' to each cluster
#'
#' @export get_ISs

get_ISs = function(x, highlight=c()) {
  highlight = get_highlight(x, highlight=highlight)
  return( sapply(highlight, get_ISs_single_cluster, x=x) )
}


get_ISs_single_cluster = function(x, cluster) {
  return(
    x %>%
      get_cov_dataframe() %>%
      dplyr::filter(labels==cluster) %>%
      dplyr::pull(IS) %>%
      unique() %>%
      length()
  )
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

get_labels = function(x, init=F) {
  py_model = get_model(x)
  tryCatch(
    expr = {
      clusters_sort = get_unique_labels(py_model, init=init)

      if (init)
        return(factor(paste("C", py_model$init_params$clusters$detach()$numpy(), sep=""), levels=clusters_sort))
      return(factor(paste("C", py_model$params$clusters$detach()$numpy(), sep=""), levels=clusters_sort))
      },
    error = function(e) return(list()) )

  check_params(x)
  if ("labels" %in% names(x$params) && !init) return(x$params$labels)

  else if ("labels_init" %in% names(x$params) && init) return(x$params$labels_init)

  cli::cli_alert_warning("No labels assignments stored in the object.")
  return(list())
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

get_unique_labels = function(x, init=FALSE) {
  if (init)
    try(expr = { labels = x$params$labels_init %>% levels(); if (!purrr::is_empty(labels)) return(labels) },
        silent = T)
  try(expr = { labels = x$params$labels %>% levels(); if (!purrr::is_empty(labels)) return(labels) },
      silent = T)

  py_model = get_model(x)
  if (init)
    tryCatch(
      expr = {
        clusters_sort = paste("C", py_model$init_params$clusters$detach()$numpy() %>% unique() %>% sort(), sep="")
        print(clusters_sort)
        return(clusters_sort) },
      error = function(e) return(0:(py_model$params$K - 1)) )

  tryCatch(
    expr = {
      clusters_sort = paste("C", py_model$params$clusters$detach()$numpy() %>% unique() %>% sort(), sep="")
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

get_unique_muts_labels = function(x, clusters=c()) {
  if (!have_vaf_df(x)) return(list())

  if (purrr::is_empty(clusters)) return(get_all_unique_muts_labels(x))

  return(
    x %>%
      get_vaf_dataframe() %>%
      dplyr::filter(labels %in% clusters) %>%
      dplyr::pull(labels_mut) %>%
      unique()
  )
}


get_all_unique_muts_labels = function(x) {
  if (!have_vaf_df(x)) return(list())
  return(
    x %>%
      get_vaf_dataframe() %>%
      tidyr::drop_na() %>%
      dplyr::pull(labels_mut) %>%
      unique()
    )
}


get_growth_rates = function(x) {
  if ("growth.rates" %in% names(x)) return(x$growth.rates) else return(data.frame())
}


check_params = function(x) {
  if (!"params" %in% names(x)) return(cli::format_warning("No params stored in the object."))
}

