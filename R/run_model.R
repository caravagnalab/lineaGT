#' Filters the input dataset.
#'
#' @description
#'
#' @param dataset Input dataset. It has to be of shape \code{NxT}, where \code{N} is the number of
#' observations and \code{T} is the dimensionality of the dataset.
#' @param columns List of columns to be selected from the dataset. If not specified,
#' all columns will be considered.
#' @param IS_values List of Insertion Sites identifier. It must be a list of length \code{N}, hence
#' for each observation a value is needed. If not specified, a general identifier \code{IS.1,...IS.N}
#' will be generated.
#' @param min_cov add
#' @param min_ccf add
#' @param k_interval add
#' @param metric add
#' @param random_state add
#' @return a dataset of the same shape as the input one, with filtered observations.
#'
#' @export filter_dataset

filter_dataset = function(dataset, columns=list(), IS=list(), min_cov=50, min_ccf=0.05, k_interval=c(5,30),
                          metric="calinski_harabasz_score", random_state=25) {
  py_pkg = reticulate::import("pylineaGT")
  x = initialize_object(K=as.integer(1), dataset=dataset, lineages=list(), columns=columns, IS_values=IS)
  x$py_model$filter_dataset(min_cov=as.integer(min_cov), min_ccf=as.numeric(min_ccf), metric=metric,
                              k_interval=as.integer(k_interval), random_state=as.integer(random_state))
  return(get_python_dataframe(x$py_model))
}


initialize_object = function(K, dataset, lineages, columns=list(), IS_values=list()) {
  py_pkg = reticulate::import("pylineaGT")
  if (purrr::is_empty(columns)) columns = dataset %>% colnames()
  if (purrr::is_empty(IS_values)) IS_values = paste("IS", 1:nrow(dataset), sep=".")
  py_model = py_pkg$mvnmm$MVNMixtureModel(K=as.integer(K), data=dataset %>% dplyr::select(all_of(columns)),
                                          lineages=lineages, IS=IS_values, columns=columns)
  return(get_object(py_model))
}


run_inference = function(x, steps=500, covariance="diag", lr=0.005, random_state=25) {
  x$py_model$fit(steps=as.integer(steps), cov_type=covariance, lr=as.numeric(lr),
                 random_state=as.integer(random_state), convergence=TRUE)
  return(update_params(x))
}


classifier = function(x) {
  x$py_model$classifier()
  return(get_object(x$py_model))
}


