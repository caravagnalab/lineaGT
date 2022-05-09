#' Retrieve the object dataframe.
#'
#' @description
#'
#' @param x a mvnmm object.
#' @return the dataset used to fit the model.
#'
#' @export get_dataframe

get_dataframe = function(x) {
  try(expr = { dataframe = x$cov.dataframe; if (!purrr::is_empty(dataframe)) return(dataframe) }, silent = T)

  py_model = get_model(x)
  return(get_python_dataframe(py_model))
}


#' Retrieve the object VAF dataframe.
#'
#' @description
#'
#' @param x a mvnmm object.
#' @return the VAF dataset.
#'
#' @export get_vaf_dataframe

get_vaf_dataframe = function(x) {
  return(x$vaf.dataframe)
}


#' Retrieve the unique labels of mutation clusters, given the coverage labels.
#'
#' @description
#'
#' @param x a mvnmm object.
#' @param clusters a vector-like variable, with the coverage labels.
#' @return vector of mutations labels.
#'
#' @export get_viber_clusters

get_viber_clusters = function(x, clusters) {
  if (purrr::is_empty(clusters)) return(get_unique_viber_labels(x))
  vaf = get_vaf_dataframe(x) %>% filter(labels %in% clusters)
  return(vaf$labels_mut %>% unique())
}


get_model = function(x) {
  if ("mvnmm" %in% class(x)) return(x$py_model)
  if ("pylineaGT.mvnmm.MVNMixtureModel" %in% class(x)) return(x)
}

