#' Retrieve the coverage dataframe.
#'
#' @description Function to retrieve the coverage dataframe used to initialize the object.
#'
#' @param x a mvnmm object.
#' @return the coverage dataset used to fit the model.
#'
#' @examples
#' if (FALSE) get_cov_dataframe(x)
#'
#' @export get_cov_dataframe

get_cov_dataframe = function(x) {
  try(expr = {
    dataframe = x$cov.dataframe; if (!purrr::is_empty(dataframe)) return(dataframe)
    },
    silent = T)

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
#' @examples
#' if (FALSE) get_vaf_dataframe(x)
#'
#' @export get_vaf_dataframe
#'

get_vaf_dataframe = function(x, label="") {
  if (paste("vaf.dataframe", label, sep=".") %in% names(x))
    return(x[[paste("vaf.dataframe", label, sep=".")]])
  return(x$vaf.dataframe)
}


#' Extract the model dimensions.
#'
#' @description Returns a vector with the dimensions of the model.
#'
#' @param x a mvnmm object.
#' @return vector of model dimensions.
#'
#' @examples
#' if (FALSE) get_dimensions(x)
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
#' if (FALSE) get_lineages(x)
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
#' if (FALSE) get_timepoints(x)
#'
#' @export get_timepoints

get_timepoints = function(x) {
  return(x$timepoints)
}


get_tp_to_int = function(x) {
  if ("tp_to_int" %in% names(x)) return(x$tp_to_int)
  return(list())
}


get_muts_fit = function(x, label="") {
  if (label == "") return(x$viber_run)
  return(x[[paste("viber_run", label, sep=".")]])
}


get_model = function(x) {
  if ("mvnmm" %in% class(x)) return(x$py_model)
  if ("pylineaGT.mvnmm.MVNMixtureModel" %in% class(x)) return(x)
}


get_color_palette = function(x, label="") {
  if (label=="") return(x$color_palette)

  if (paste("color_palette", label, sep=".") %in% names(x))
    return(x[[paste("color_palette", label, sep=".")]])
  else
    return(x$color_palette)
}


get_trees = function(x, label="") {
  if (!any(grepl("trees", x %>% names)))
    return(list())
  if (label=="")
    return(x$trees)
  return(x[[paste("trees", label, sep=".")]])
}


get_adj = function(tree.k) {
  return(
    tree.k$adj_mat
  )
}
