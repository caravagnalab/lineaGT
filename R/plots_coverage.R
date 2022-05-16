#' 2D scatterplot and density
#'
#' @description Function to plot the scatterplots of the coverage, one timepoint against the other, together with the Gaussian density.
#' @param x a mvnmm object.
#' @param plot_density a Boolean. If set to FALSE, the Gaussian density will not be displayed.
#' @param facet a Boolean. If set to TRUE, the plot will be faceted on the labels.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#'
#' @return a list of plots, one for each timepoints combination.
#'
#' @examples
#' if (FALSE) plots = plot_scatter_density(x)
#'
#' @import ggplot2
#' @importFrom purrr is_empty map
#' @importFrom dplyr select all_of
#' @export plot_scatter_density


plot_scatter_density = function(x, plot_density=T, facet=FALSE, highlight=c()) {
  dataset = x %>% get_cov_dataframe() %>% long_to_wide_cov()

  color_palette = highlight_palette(x$color_palette, highlight)
  if (purrr::is_empty(highlight)) highlight = get_unique_labels(x)

  if (plot_density) density = compute_density(x) else density = NULL

  combinations = get_pairs(dataset, columns=x %>% get_dimensions())
  max_val = max(dataset %>% dplyr::select(dplyr::all_of(x %>% get_dimensions()))) +10

  p = list()
  for (t1_t2 in combinations$pair_name) {
    xy = strsplit(t1_t2, ":")[[1]]
    p[[t1_t2]] = plot_2D(x, xy[1], xy[2], color_palette, highlight, dens=density,
                         facet=facet, ylim=c(0,max_val), xlim=c(0,max_val))
  }
  return(p)
}


plot_2D = function(x, dim1, dim2, color_palette, highlight, dens=NULL, facet=F, ...) {
  inputs = eval(substitute(alist(...))) %>% purrr::map(as.list)

  pl = ggplot2::ggplot() +
    geom_point(data=x %>% get_cov_dataframe() %>% long_to_wide_cov(),
               aes_string(x=dim1, y=dim2, color="labels"),
               alpha=.4, size=.8) +
    scale_color_manual(values=color_palette, breaks=highlight) +
    labs(color="Clusters") +
    xlab(split_to_camelcase(dim1)) +
    ylab(split_to_camelcase(dim2)) +
    my_ggplot_theme() +
    theme(legend.position="bottom") +
    coord_fixed(ratio=1)

  try({ pl = pl + ylim(inputs$ylim) }, silent=T)
  try({ pl = pl + xlim(inputs$xlim) }, silent=T)

  if (!is.null(dens)) pl = pl +
    geom_density_2d(data=dens %>% filter(labels %in% highlight),
                    aes_string(x=dim1, y=dim2, color="labels"),
                    inherit.aes=F, contour_var="ndensity", size=.1)

  if (facet) pl = pl + facet_wrap(~labels, nrow=1) + theme(legend.position="none")
  return(pl)
}


#' Histogram of the marginal distribution of each dimension
#'
#' @description Function to plot the marginal distribution of the coverage, for each timepoint, colored by cluster.
#' @param x a mvnmm object.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param binwidth numeric value representing the histogram binwidth.
#'
#' @import ggplot2
#'
#' @examples
#' if (FALSE) plot_marginal(x)
#'
#' @export plot_marginal


plot_marginal = function(x, highlight=c(), binwidth=5) {
  color_palette = highlight_palette(x$color_palette, highlight)
  if (purrr::is_empty(highlight)) highlight = get_unique_labels(x)

  dd = x %>% get_cov_dataframe() %>% filter(labels %in% highlight)

  p = dd %>% ggplot() +
    geom_histogram(aes(x=coverage, fill=labels), position="identity", alpha=.7, binwidth=binwidth) +
    scale_fill_manual(values=color_palette, breaks=highlight) +
    facet_grid(timepoints ~ labels) +
    ylab("Counts") +
    xlab("Coverage") +
    labs(fill="Clusters") +
    my_ggplot_theme()

  return(p)
}
