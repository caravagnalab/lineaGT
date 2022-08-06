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
#' @importFrom dplyr select all_of group_by filter pull
#' @importFrom purrr is_empty
#'
#' @export plot_scatter_density

plot_scatter_density = function(x, plot_density=T, highlight=c(), min_frac=0) {

  dataset = x %>% get_cov_dataframe() %>% long_to_wide_cov()

  highlight = get_highlight(x, min_frac=min_frac, highlight=highlight)
  color_palette = highlight_palette(x, highlight)

  if (plot_density) density = compute_density(x) else density = NULL

  combinations = get_pairs(dataset, columns=x %>% get_dimensions())
  max_val = max(dataset %>% dplyr::select(dplyr::all_of(x %>% get_dimensions()))) + 10

  p = list()
  for (t1_t2 in combinations$pair_name) {
    xy = strsplit(t1_t2, ":")[[1]]
    p[[t1_t2]] = plot_2D(x,
                         xy[1],
                         xy[2],
                         color_palette,
                         highlight,
                         dens=density,
                         ylim=c(0,max_val),
                         xlim=c(0,max_val))
  }
  return(p)
}


plot_2D = function(x, dim1, dim2, color_palette, highlight, dens=NULL, ...) {
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
#' @importFrom purrr is_empty
#' @importFrom dplyr mutate filter inner_join group_by ungroup
#' @importFrom tidyr unnest
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' if (FALSE) plot_marginal(x)
#'
#' @export plot_marginal

plot_marginal = function(x, min_frac=0, highlight=c(), binwidth=10, show_dens=T, timepoints_to_int=list()) {

  if (purrr::is_empty(timepoints_to_int)) timepoints_to_int = map_timepoints_int(x)

  tp = timepoints_to_int %>% unlist() %>% sort() %>% names()

  highlight = get_highlight(x, min_frac, highlight)
  color_palette = highlight_palette(x, highlight)

  dd = x %>%
    get_cov_dataframe() %>%
    dplyr::mutate(timepoints=factor(timepoints, levels=tp)) %>%
    dplyr::filter(labels %in% highlight)

  means = x %>%
    get_mean_long() %>%
    dplyr::mutate(timepoints=factor(timepoints, levels=tp)) %>%
    dplyr::filter(labels %in% highlight)

  vars = x %>%
    get_variance_long() %>%
    dplyr::mutate(timepoints=factor(timepoints, levels=tp)) %>%
    dplyr::filter(labels %in% highlight)

  weights = (get_ISs(x) / sum(get_ISs(x))) %>%
    data.frame() %>% setNames("weights") %>%
    tibble::rownames_to_column(var="labels")

  params = dplyr::inner_join(means, vars, by=c("labels", "timepoints", "lineage")) %>%
    dplyr::inner_join(weights, by="labels")

  p = list()
  lineages = x %>% get_lineages()
  for (ll in lineages) {
    dens.df = params %>%
      dplyr::filter(lineage==ll) %>%
      dplyr::group_by(timepoints, labels) %>%
      dplyr::mutate(dens=list(rnorm(1000, mean=mean_cov, sd=sigma))) %>%
      tidyr::unnest(dens) %>%
      dplyr::filter(dens>=0) %>%
      dplyr::ungroup()

    p[[ll]] = dd %>%
      dplyr::filter(lineage==ll) %>%
      ggplot() +
      geom_histogram(aes(x=coverage, fill=labels, y=..count../sum(..count..)), position="identity", alpha=.7, binwidth=binwidth, inherit.aes=F) +
      scale_fill_manual(values=color_palette, breaks=highlight) +
      scale_color_manual(values=color_palette, breaks=highlight) +
      facet_grid(timepoints ~ labels) +
      ylab("Counts") +
      xlab("Coverage") +
      labs(fill="Clusters", subtitle=ll) +
      my_ggplot_theme() + guides(color="none")

    if (show_dens)
      p[[ll]] = p[[ll]] + geom_density(data=dens.df, aes(x=dens, y=..count../sum(..count..)*binwidth, color=labels),
                                       # color="#424949",
                                       linetype="solid", size=.1, inherit.aes=F)
  }

  return(p)
}
