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
    geom_point(
      data=long_to_wide_cov(get_cov_dataframe(x)),
      # data = dens %>% filter(labels %in% highlight),
      aes_string(x=dim1, y=dim2, color="labels"),
      alpha=.8, size=.8) +
    scale_color_manual(values=color_palette, breaks=highlight) +
    labs(color="Clusters") +
    xlab(split_to_camelcase(dim1)) +
    ylab(split_to_camelcase(dim2)) +
    my_ggplot_theme() +
    theme(legend.position="bottom")
    # coord_fixed(ratio=1)

  try({ pl = pl + ylim(inputs$ylim) }, silent=T)
  try({ pl = pl + xlim(inputs$xlim) }, silent=T)

  if (!is.null(dens))
    pl = pl +
      stat_density_2d(data=dens %>% filter(labels %in% highlight),
                      aes_string(x=dim1, y=dim2, alpha="..level..", fill="labels"),
                      geom="polygon", inherit.aes=F, contour_var="ndensity", bins=10) +

      scale_fill_manual(values=color_palette, breaks=highlight) +
      scale_alpha_continuous(range=c(0.01,0.3)) +
      labs(fill="Clusters") + guides(alpha="none", color="none")

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

plot_marginal = function(x,
                         min_frac=0,
                         highlight=c(),
                         binwidth=10,
                         show_dens=T,
                         timepoints_to_int=list(),
                         single_plot=F) {

  timepoints_to_int = map_timepoints_int(x)

  tp = timepoints_to_int %>% unlist() %>% sort() %>% names()

  highlight = get_highlight(x, min_frac, highlight)
  color_palette = highlight_palette(x, highlight)

  dd = x %>%
    get_cov_dataframe() %>%
    dplyr::mutate(timepoints=factor(timepoints, levels=tp))

  if (show_dens) {
    means = x %>%
      get_mean_long() %>%
      dplyr::mutate(timepoints=factor(timepoints, levels=tp))

    vars = x %>%
      get_variance_long() %>%
      dplyr::mutate(timepoints=factor(timepoints, levels=tp))

    weights = (get_ISs(x) / sum(get_ISs(x))) %>%
      data.frame() %>% setNames("weights") %>%
      tibble::rownames_to_column(var="labels")

    params = dplyr::inner_join(means, vars, by=c("labels", "timepoints", "lineage")) %>%
      dplyr::inner_join(weights, by="labels") %>%
      dplyr::mutate(labels=factor(labels, levels=get_unique_labels(x)))
  }


  if (single_plot) {
    max_val = max(dd$coverage)
    dens.df = params %>%
      dplyr::group_by(lineage, timepoints, labels) %>%
      dplyr::mutate(dens=list(weights * dnorm(1:max_val, mean=mean_cov, sd=sigma)),
                    coverage=list(1:max_val)) %>%
      tidyr::unnest(c(dens, coverage)) %>%
      dplyr::filter(dens>=0) %>%
      dplyr::ungroup()

    dens.all = dens.df %>%
      group_by(timepoints, lineage, coverage) %>%
      dplyr::summarise(dens=sum(dens)) %>%
      unique()

    p = ggplot() +
      geom_density(data=dens.all,
                   aes(x=coverage, y=dens*binwidth),
                   stat="identity", inherit.aes=F, fill="#FFFFFF00") +
      geom_density(data=dens.df,
                   aes(x=coverage, y=dens*binwidth, fill=labels), alpha=.3,
                   stat="identity", inherit.aes=F, color="#FFFFFF00") +
      geom_rug(data=filter(dd, labels %in% highlight), aes(x=coverage, color=labels)) +

      scale_fill_manual(values=color_palette, breaks=highlight) +
      scale_color_manual(values=color_palette, breaks=highlight) +
      ylab("Counts") + xlab("Coverage") +
      labs(fill="Clone") +
      my_ggplot_theme() +
      ylab("Density") + guides(color="none") +
      guides(fill=guide_legend(override.aes=list(alpha=1)))

    if ((dd$lineage %>% unique() %>% length()) >= 1)
      p = p + facet_grid(timepoints ~ lineage, scale="free_y") else
        p = p + facet_grid(timepoints ~ ., scale="free_y")

    return(p + xlim(0, max(x$cov.dataframe$coverage)))
  }


  p = list()
  lineages = x %>% get_lineages()
  for (ll in lineages) {

    p[[ll]] = dd %>%
      dplyr::filter(labels %in% highlight) %>%
      dplyr::filter(lineage==ll) %>%
      ggplot() +
      geom_histogram(aes(x=coverage,
                         y=..count../sum(..count..)/binwidth,
                         fill=labels),
                     position="identity", alpha=1, binwidth=binwidth, color="#FFFFFF00") +
      scale_fill_manual(values=color_palette, breaks=highlight) +
      facet_grid(timepoints ~ labels) +
      ylab("Counts") + xlab("Coverage") +
      labs(fill="Clusters", subtitle=ll) +
      my_ggplot_theme()

    if (show_dens) {
      dens.df = params %>%
        dplyr::filter(lineage==ll) %>%
        dplyr::group_by(timepoints, labels) %>%
        # dplyr::mutate(dens=list(rnorm(1000, mean=mean_cov, sd=sigma))) %>%
        dplyr::mutate(dens=list(dnorm(1:max(dd$coverage), mean=mean_cov, sd=sigma)),
                      coverage=list(1:max(dd$coverage))) %>%
        tidyr::unnest(c(coverage, dens)) %>%
        # dplyr::filter(dens>=0) %>%
        dplyr::ungroup()

      p[[ll]] = p[[ll]] +
        geom_density(data=dens.df,
                     aes(x=coverage, y=dens, fill=labels),
                     stat="identity",
                     position="identity", color="#FFFFFF00", inherit.aes=F, alpha=.3) +
        ylab("Density") + guides(color="none")
    }

  }

  return(p)
}
