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

get_cov_dataframe = function(x, verbose=T) {
  try(expr = {
      dataframe = x$cov.dataframe
      if (!purrr::is_empty(dataframe)) return(dataframe) else if (verbose)
        return(cli::cli_alert_warning("No coverage dataframe loaded.")) else
        return(NULL)
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

get_vaf_dataframe = function(x, verbose=T) {
  if ("vaf.dataframe" %in% names(x)) return(x$vaf.dataframe)

  if (verbose) cli::cli_alert_warning("No VAF dataframe loaded.")
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

get_dimensions = function(x, verbose=T) {
  if ("dimensions" %in% names(x)) return(x$dimensions)

  if (verbose) cli::cli_alert_warning("No dimensions values.")
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

get_lineages = function(x, verbose=T) {
  if ("lineages" %in% names(x)) return(x$lineages)

  if (verbose) cli::cli_alert_warning("No lineages values.")
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


get_ISs_list = function(x, cls=c(), mutations=F) {
  if (length(cls)==0) cls = x %>% get_highlight(mutations=mutations)
  cls.tmp = list("clonal"=get_highlight(x, mutations = F), "subcl"=get_all_unique_muts_labels(x))

  clonal = lapply(cls.tmp$clonal, function(cc)
    x %>% get_cov_dataframe() %>%
      dplyr::filter(labels==cc) %>%
      dplyr::pull(IS) %>%
      unique()) %>% setNames(cls.tmp$clonal)

  subcl = lapply(cls.tmp$subcl, function(ss)
    x %>% get_vaf_dataframe() %>%
      dplyr::filter(labels_mut==ss) %>%
      dplyr::pull(IS) %>%
      unique()) %>% setNames(cls.tmp$subcl)

  list_ISs = c(clonal, subcl)

  return(list_ISs[cls])
}


get_tp_to_int = function(x, verbose=T) {
  if ("tp.to.int" %in% names(x)) return(x$tp.to.int)

  if (verbose) cli::cli_alert_warning("Timepoints not mapped to integer values.")
  return(list())
}


get_muts_fit = function(x, verbose=T) {
  if ("x.muts" %in% names(x)) return(x$x.muts)

  if (verbose) cli::cli_alert_warning("No fitted binomial clusters object.")
  return(list())
}


get_model = function(x) {
  if ("mvnmm" %in% class(x)) return(x$py_model)
  if ("pylineaGT.mvnmm.MVNMixtureModel" %in% class(x)) return(x)
}


get_color_palette = function(x, verbose=T) {
  if ("color.palette" %in% names(x)) return(x$color.palette)

  if (verbose) cli::cli_alert_warning("No color palette.")
  return(list())
}


get_trees = function(x, verbose=T) {
  if ("x.trees" %in% names(x)) return(x$x.trees)

  if (verbose) cli::cli_alert_warning("No fitted phylogenies object.")
  return(list())
}


get_adj = function(tree.k) {
  return(
    tree.k$adj_mat
  )
}


get_pop_df = function(x, verbose=T) {
  if (have_pop_df(x))
    return(x$population.df)

  if (verbose) cli::cli_alert_warning("No stored population dataframe!")
  return(NULL)
}
