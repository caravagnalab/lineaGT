#' Filters the input dataset.
#'
#' @description Function used to filter observations, i.e. ISs, in the input dataframe with coverage values.
#'
#' @param cov.df Input coverage dataset. It must have at least the columns \code{coverage}, \code{timepoints},
#' \code{lineage}, \code{IS}, with the coverage values, timepoint, lineage and IS, respectively.
#' @param min_cov add
#' @param min_frac add
#' @param k_interval add
#' @param metric add
#' @param seed add
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
                          k_interval=c(10,20),
                          metric="calinski_harabasz_score",
                          seed=5) {

  cli::cli_process_start(paste0("Filtering the input dataset with minimum coverage ", min_cov, " and minimum clusters fraction ", min_frac,"."))

  # if (!is.factor(cov.df$timepoints))
  #   cov.df = cov.df %>%
  #     dplyr::mutate(timepoints=as.character(timepoints))

  cov.df = cov.df %>%
    check_cov_dimensions() %>%
    dplyr::group_by(IS) %>%
    dplyr::filter(any(coverage>=min_cov)) %>%
    dplyr::ungroup()

  if(nrow(cov.df)==0) {
    cli::cli_process_failed(msg_failed="The filtered dataset contains 0 ISs.")
    return(cov.df)
  }

  if (min_frac == 0) {
    cli::cli_process_done()
    return(cov.df)
  }

  max_k = cov.df %>% check_max_k()
  k_interval = check_k_interval(k_interval, max_k)

  py_pkg = reticulate::import("pylineaGT")

  py_model = initialize_object(K=as.integer(1), cov.df=cov.df, py_pkg=py_pkg, return_model=TRUE)

  py_model$filter_dataset(min_cov=as.integer(0),
                          min_frac=as.numeric(min_frac),
                          metric=metric,
                          k_interval=as.integer(k_interval),
                          seed=as.integer(seed))

  cli::cli_process_done()
  return(
      get_python_dataframe(py_model)
    )
}

