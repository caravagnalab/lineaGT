#' Filters the input dataset.
#'
#' @description
#'
#' @param cov.df Input coverage dataset. It must have at least the columns \code{coverage}, \code{timepoints},
#' \code{lineage}, \code{IS}, with the coverage values, timepoint, lineage and IS, respectively.
#' @param min_cov add
#' @param min_frac add
#' @param k_interval add
#' @param metric add
#' @param random_state add
#' @return a dataset of the same shape as the input one, with filtered observations.
#'
#' @importFrom reticulate import
#' @importFrom magrittr %>%
#' @importFrom dplyr select starts_with
#'
#' @export filter_dataset

filter_dataset = function(cov.df,
                          min_cov=5,
                          min_frac=0.05,
                          k_interval=c(5,30),
                          metric="calinski_harabasz_score",
                          random_state=25) {

  py_pkg = reticulate::import("pylineaGT")
  x = initialize_object(K=as.integer(1), cov.df=cov.df, py_pkg)
  x$py_model$filter_dataset(min_cov=as.integer(min_cov),
                            min_ccf=as.numeric(min_frac),
                            metric=metric,
                            k_interval=as.integer(k_interval),
                            random_state=as.integer(random_state))
  return(get_python_dataframe(x$py_model))
}

