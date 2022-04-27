#' Plot scatter and density
#'
#' @description Function to plot the scatterplots of the coverage, one timepoint against the other, together with the Gaussian density.
#' @param obj a mvnmm object.
#' @param plot_density a Boolean. If set to FALSE, the Gaussian density will not be displayed.
#' @param facet a Boolean. If set to TRUE, the plot will be faceted on the labels.
#'
#' @return a list of plots, one for each timepoints combination.
#'
#' @import ggplot2
#' @importFrom purrr is_empty map
#' @importFrom dplyr select all_of
#' @export plot_scatter_density
#'
#' @examples
#' plot_list = plot_scatter_density(obj, plot_density=TRUE, facet=FALSE)
#' plot_list$`cov_early_B:cov_early_MNC`

plot_scatter_density = function(obj, plot_density=T, facet=FALSE, highlight=c()) {
  py_model = obj$py_model
  timepoints = obj$dimensions
  params = obj$params

  dataset = obj$dataframe

  color_palette = highlight_palette(obj$color_palette, highlight)
  if (purrr::is_empty(highlight)) highlight = get_unique_labels(obj)

  if (plot_density) density = lineaGT::compute_density(obj) else density = NULL

  combinations = get_pairs(dataset, columns=obj$dimensions)
  max_val = max(dataset %>% dplyr::select(dplyr::all_of(obj$dimensions))) +10

  p = list()
  for (t1_t2 in combinations$pair_name) {
    xy = strsplit(t1_t2, ":")[[1]]
    p[[t1_t2]] = plot_2D(obj, xy[1], xy[2], color_palette, highlight, dens=density,
                         facet=facet, ylim=c(0,max_val), xlim=c(0,max_val))
  }
  return(p)
}


plot_2D = function(obj, dim1, dim2, color_palette, highlight, dens=NULL, facet=F, ...) {
  inputs = eval(substitute(alist(...))) %>% purrr::map(as.list)

  pl = ggplot2::ggplot() +
    geom_point(data=obj$dataframe, aes_string(x=dim1, y=dim2, color="labels"), alpha=.4, size=.8) +
    scale_color_manual(values=color_palette) + labs(color="Clusters") +
    xlab(split_to_camelcase(dim1)) + ylab(split_to_camelcase(dim2)) +
    my_ggplot_theme() + theme(legend.position="bottom") + coord_fixed(ratio=1)

  try({ pl = pl + ylim(inputs$ylim) }, silent=T)
  try({ pl = pl + xlim(inputs$xlim) }, silent=T)

  if (!is.null(dens)) pl = pl + geom_density_2d(data=dens %>% filter(labels %in% highlight),
                                              aes_string(x=dim1, y=dim2, color="labels"),
                                              inherit.aes=F, contour_var="ndensity", size=.1) +
      scale_color_manual(values=obj$color_palette)

  if (facet) pl = pl + facet_wrap(~labels, nrow=1) + theme(legend.position="none")
  return(pl)
}

