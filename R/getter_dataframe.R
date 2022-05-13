#' Retrieve the coverage dataframe.
#'
#' @description Function to retrieve the coverage dataframe used to initialize the object.
#'
#' @param x a mvnmm object.
#' @return the coverage dataset used to fit the model.
#'
#' @export get_cov_dataframe

get_cov_dataframe = function(x) {
  try(expr = { dataframe = x$cov.dataframe; if (!purrr::is_empty(dataframe)) return(dataframe) }, silent = T)

  py_model = get_model(x)
  return(get_python_dataframe(py_model))
}


#' Retrieve the mutations dataframe.
#'
#' @description Function to retrieve the mutations dataframe used to initialize the object.
#'
#' @param x a mvnmm object.
#' @return the mutations dataset.
#'
#' @export get_vaf_dataframe
#'
get_vaf_dataframe = function(x) {
  return(x$vaf.dataframe)
}


get_model = function(x) {
  if ("mvnmm" %in% class(x)) return(x$py_model)
  if ("pylineaGT.mvnmm.MVNMixtureModel" %in% class(x)) return(x)
}


#' Extract the model dimensions.
#'
#' @description Returns a vector with the dimensions of the model.
#'
#' @param x a mvnmm object.
#' @return vector of model dimensions.
#'
#' @examples
#' get_dimensions(x)
#'
#' @export get_dimensions

get_dimensions = function(x) {
  return(x$dimensions)
}

#' Extract the data lineages.
#'
#' @description Returns a vector with the lineages of the input data.
#'
#' @param x a mvnmm object.
#' @return vector of data lineages.
#'
#' @examples
#' get_lineages(x)
#'
#' @export get_lineages

get_lineages = function(x) {
  return(x$lineages)
}


#' Extract the data timepoints.
#'
#' @description Returns a vector with the timepoints of the input data.
#'
#' @param x a mvnmm object.
#' @return vector of data timepoints.
#'
#' @examples
#' get_timepoints(x)
#'
#' @export get_timepoints

get_timepoints = function(x) {
  return(x$timepoints)
}
