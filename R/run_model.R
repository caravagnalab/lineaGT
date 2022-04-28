#' Creates an object of class mvnmm.
#'
#' @description Add
#'
#' @param K Input maximum number of clusters
#' @param dataset Input dataset. It has to be of shape NxT, where N is the number of
#' observations and T is the dimensionality of the dataset.
#' @param columns List of columns to be selected from the dataset. If not specified,
#' all columns will be considered.
#' @param IS_values List of Insertion Sites identifier. It must be a list of length N, hence
#' for each observation a value is needed. If not specified, a general identifier IS.1,...IS.N
#' will be generated.
#' @return a LineaGT object, containing the input dataset, annotated with IS_values,
#' N, K, T specific of the dataset, the input IS and column names, a list params that
#' will contain the inferred parameters, the python object
#'
#' @importFrom dplyr filter mutate select group_by inner_join rename_with case_when all_of ungroup
#' @importFrom tidyr separate unite pivot_wider pivot_longer tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace_all
#' @importFrom MASS mvrnorm
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom reshape2 melt dcast
#' @importFrom grDevices gray
#' @importFrom stats coef lm runif
#' @importFrom utils combn read.csv
#' @importFrom VIBER variational_fit choose_clusters
#'
#' @export mixture_model
#' @export filter_dataset
#' @export run_inference
#' @export classifier


mixture_model = function(K, dataset, lineages, columns=list(), IS_values=list()) {
  py_pkg = reticulate::import("pylineaGT")

  if (purrr::is_empty(columns)) columns = dataset %>% colnames()
  if (purrr::is_empty(IS_values)) IS_values = paste("IS", 1:nrow(dataset), sep=".")
  py_model = py_pkg$mvnmm$MVNMixtureModel(K=as.integer(K), lineages=lineages,
                                          data=dataset %>% dplyr::select(all_of(columns)),
                                          IS=IS_values, columns=columns)
  return(get_object(py_model))
}


filter_dataset = function(dataset, columns, IS, min_cov=50, min_ccf=0.05, k_interval=c(5,30),
                          metric="calinski_harabasz_score", random_state=25) {
  py_pkg = reticulate::import("pylineaGT")
  py_model = py_pkg$mvnmm$MVNMixtureModel(K=as.integer(1), lineages=list(),
                                          data=dataset %>% dplyr::select(columns),
                                          columns=columns, IS=IS)
  py_model$filter_dataset(min_cov=as.integer(min_cov),
                          min_ccf=as.numeric(min_ccf),
                          k_interval=as.integer(k_interval),
                          metric=metric,
                          random_state=as.integer(random_state))
  return(get_dataset(py_model))
}


run_inference = function(obj, steps=500, covariance="diag", lr=0.005, random_state=25) {
  obj$py_model$fit(steps=as.integer(steps),
                   cov_type=covariance,
                   lr=as.numeric(lr),
                   random_state=as.integer(random_state))
  return(update_params(obj))
}


classifier = function(obj) {
  obj$py_model$classifier()
  return(get_object(obj$py_model))
}


