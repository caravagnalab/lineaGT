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
#' @param min_frac add
#' @param k_interval add
#' @param metric add
#' @param random_state add
#' @return a dataset of the same shape as the input one, with filtered observations.
#'
#' @importFrom reticulate import
#'
#' @export filter_dataset

filter_dataset = function(dataset,
                          min_cov=50,
                          min_frac=0.05,
                          k_interval=c(5,30),
                          metric="calinski_harabasz_score",
                          random_state=25) {

  py_pkg = reticulate::import("pylineaGT")
  x = initialize_object(K=as.integer(1), dataset=dataset)
  x$py_model$filter_dataset(min_cov=as.integer(min_cov),
                            min_ccf=as.numeric(min_frac),
                            metric=metric,
                            k_interval=as.integer(k_interval),
                            random_state=as.integer(random_state))
  return(get_python_dataframe(x$py_model))
}

