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
#' @import ggplot2
#' @importFrom purrr is_empty map
#' @importFrom dplyr select all_of
#' @export plot_scatter_density


plot_scatter_density = function(x, plot_density=T, facet=FALSE, highlight=c()) {
  py_model = x$py_model
  timepoints = x$dimensions
  params = x$params

  dataset = x$dataframe

  color_palette = highlight_palette(x$color_palette, highlight)
  if (purrr::is_empty(highlight)) highlight = get_unique_labels(x)

  if (plot_density) density = compute_density(x) else density = NULL

  combinations = get_pairs(dataset, columns=x$dimensions)
  max_val = max(dataset %>% dplyr::select(dplyr::all_of(x$dimensions))) +10

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
    geom_point(data=x$dataframe, aes_string(x=dim1, y=dim2, color="labels"), alpha=.4, size=.8) +
    scale_color_manual(values=color_palette) +
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
#' @export plot_marginal


plot_marginal = function(x, highlight=c(), binwidth=5) {
  weights = get_weights(x)
  max_value = get_dataframe(x) %>% dplyr::select(starts_with("cov")) %>% max()
  data = x %>% get_dataframe()

  color_palette = highlight_palette(x$color_palette, highlight)
  if (purrr::is_empty(highlight)) highlight = get_unique_labels(x)

  dd = x %>% get_dataframe() %>%
    dplyr::select(starts_with("cov"), "IS", starts_with("labels")) %>%
    tidyr::pivot_longer(cols=dplyr::starts_with("cov"), names_to="timepoints_lineage", values_to="cov") %>%
    tidyr::separate(timepoints_lineage, into=c("else","timepoints","lineage"), sep="\\_|\\.") %>%
    mutate("else"=NULL, timepoints=split_to_camelcase(timepoints)) %>%
    filter(labels %in% highlight)

  p = dd %>% ggplot() +
    geom_histogram(aes(x=cov, fill=labels), position="identity", alpha=.7, binwidth=binwidth) +
    scale_fill_manual(values=color_palette[highlight]) +
    facet_grid(timepoints ~ labels) +
    ylab("Counts") +
    xlab("Coverage") +
    labs(fill="Clusters") +
    my_ggplot_theme()

  return(p)
}
