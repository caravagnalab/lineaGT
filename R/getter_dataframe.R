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
      dataframe = x$cov.dataframe
      if (!purrr::is_empty(dataframe)) return(dataframe) else return(cli::cli_alert_warning("No coverage dataframe loaded."))
    }, silent = T)

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

get_vaf_dataframe = function(x) {
  if ("vaf.dataframe" %in% names(x)) return(x$vaf.dataframe)

  cli::cli_alert_warning("No VAF dataframe loaded.")
  return(list())
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
  if ("dimensions" %in% names(x)) return(x$dimensions)

  cli::cli_alert_warning("No dimensions values.")
  return(list())
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
  if ("lineages" %in% names(x)) return(x$lineages)

  cli::cli_alert_warning("No lineages values.")
  return(list())
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
  if ("timepoints" %in% names(x)) return(x$timepoints)

  cli::cli_alert_warning("No timepoints values.")
  return(list())
}


get_tp_to_int = function(x) {
  if ("tp.to.int" %in% names(x)) return(x$tp.to.int)

  cli::cli_alert_warning("Timepoints not mapped to integer values.")
  return(list())
}


get_muts_fit = function(x) {
  if ("x.muts" %in% names(x)) return(x$x.muts)

  cli::cli_alert_warning("No fitted binomial clusters object.")
  return(list())
}


get_model = function(x) {
  if ("mvnmm" %in% class(x)) return(x$py_model)
  if ("pylineaGT.mvnmm.MVNMixtureModel" %in% class(x)) return(x)
}


get_color_palette = function(x) {
  if ("color.palette" %in% names(x)) return(x$color.palette)

  cli::cli_alert_warning("No color palette.")
  return(list())
}


get_trees = function(x) {
  if ("x.trees" %in% names(x)) return(x$x.trees)

  cli::cli_alert_warning("No fitted phylogenies object.")
  return(list())
}


get_adj = function(tree.k) {
  return(
    tree.k$adj_mat
  )
}
