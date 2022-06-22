#' Muller plot
#'
#' @description Function to visualize the mullerplot for the fitted object.
#'
#' @param x a mvnmm object.
#' @param which string among \code{"frac","pop","fitness"} determining whether to plot the coverage
#' normalized in \code{[0,1]}, as absolute clone abundance, or with each clone colored by the growth rate,
#' computed assuming a exponential growth.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param timepoints_to_int a list to map each \code{timepoint} value to an integer.
#' @param mutations Boolean. If set to \code{TRUE}, also the clusters of mutations will be visualized.
#' @param single_clone Boolean. If \code{mutations} and \code{single_clone} are set to \code{TRUE}, only the clones
#' reported in \code{highlight} and the respective subclones will be visualised.
#' @param label a character corresponding to the label of the run to visualize.#' @param legend.pos position of the legend. If set to \code{"none"}, the legend is not shown.
#' @param wrap Boolean. If set to \code{TRUE}, a single plot with the mullerplots for each lineage will be returned.
#'
#' @examples
#' if (FALSE) plot_mullerplot(x, wrap=T)
#'
#' @import ggplot2
#' @import ggmuller
#' @importFrom patchwork wrap_plots
#'
#' @export plot_mullerplot

plot_mullerplot = function(x, which="frac", highlight=c(), min_frac=0,
                           timepoints_to_int=list("init"=0,"early"=60,"mid"=140,"late"=280),
                           mutations=F, single_clone=T, label="",
                           legend.pos="right", wrap=F) {
  highlight.cov = get_highlight(x, min_frac=min_frac, highlight=highlight)
  highlight = get_highlight(x, min_frac, highlight.cov, mutations=mutations, label=label)
  color_palette = highlight_palette(x, highlight, label)

  pop_df = get_muller_pop(x, mutations=mutations, map_tp_time=timepoints_to_int, label=label)
  edges_df = get_muller_edges(x, mutations=mutations, label=label)

  if (single_clone && mutations) {
    pop_df = pop_df %>% filter_muller_df(highlight=highlight.cov)
    edges_df = edges_df %>% filter_muller_df(highlight=highlight.cov)
  }

  timepoints = x %>% get_dimensions()
  lineages = x %>% get_lineages()

  plot_list = list()
  for (ll in lineages) {
    tp = timepoints[grep(pattern=ll, x=timepoints)]
    if (length(tp) != 0) {
      pop_ll = pop_df %>% dplyr::filter(Lineage==ll)
      mullerdf_ll = ggmuller::get_Muller_df(edges_df, pop_ll)

      if (which == "frac" || which == "")
        plot_list[[ll]] = mullerplot_util(mullerdf_ll,
                                          y="Frequency",
                                          fill="Identity",
                                          lineage=ll,
                                          highlight=highlight,
                                          color_palette=color_palette,
                                          legend.pos=legend.pos)
      if (which == "pop" || which == "")
        plot_list[[ll]] = mullerplot_util(mullerdf_ll %>% ggmuller::add_empty_pop(),
                                          y="Population",
                                          highlight=highlight,
                                          fill="Identity",
                                          color_palette=color_palette,
                                          lineage=ll,
                                          legend.pos=legend.pos)
      if (which == "fitness"){
        exp_limits = c(min(pop_df$lm_r), max(pop_df$lm_r))
        plot_list[[ll]] = mullerplot_util(mullerdf_ll,
                                          y="Frequency",
                                          fill="lm_r",
                                          highlight=highlight,
                                          color_palette=color_palette,
                                          lineage=ll,
                                          legend.pos=legend.pos,
                                          exp_limits=exp_limits)
      }
    }
  }

  if (wrap) return(patchwork::wrap_plots(plot_list, guides="collect", ncol=2))
  return(plot_list)
}


mullerplot_util = function(mullerdf, y, fill, lineage, color_palette, highlight,
                           legend.pos="right", exp_limits=NULL) {
  if (fill=="Identity")
    pl = mullerdf %>% ggplot() +
      geom_area(aes_string(x="Generation", y=y, group="Group_id", fill="Identity", colour="Identity"), alpha=.9) +
      geom_vline(xintercept=mullerdf$Generation %>% unique(), linetype="dashed") +
      guides(linetype="none", color="none") +
      scale_fill_manual(name="Clusters", values=color_palette, na.value="transparent", breaks=highlight) +
      scale_color_manual(values=color_palette, na.value="transparent", breaks=highlight) +
      xlab("Time") +
      labs(title=split_to_camelcase(lineage)) +
      my_ggplot_theme(legend.pos=legend.pos)

  if (fill == "lm_r")
    pl = mullerdf %>% ggplot() +
      geom_area(aes_string(x="Generation", y=y, group="Group_id", fill="lm_r"), alpha=.9) +
      geom_vline(xintercept=mullerdf$Generation %>% unique(), linetype="dashed") +
      guides(linetype="none", color="none") +
      xlab("Time") + labs(title=split_to_camelcase(lineage), fill="Exp rate") +
      my_ggplot_theme(legend.pos=legend.pos) +
      scale_fill_gradient2(mid="white", low="blue", high="red", limits=exp_limits)

  return(pl)
}



#' Exponential fitting
#'
#' @description Function to visualize the growht of each clone and lineage.
#'
#' @param x a mvnmm object.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param facet Boolean. If set to \code{TRUE}, the plot will be faceted against the clusters.
#' @param mutations Boolean. If set to \code{TRUE}, the growth will be visualize for each cluster of mutations.
#' @param label a character corresponding to the label of the run to visualize.
#' @param timepoints_to_int a list to map each \code{timepoint} value to an integer.
#'
#' @examples
#' if (FALSE) plot_exp_fit(x)
#'
#' @import ggplot2
#'
#' @export plot_exp_fit

plot_exp_fit = function(x, highlight=c(), min_frac=0, facet=F, mutations=F, label="",
                        timepoints_to_int=list("init"=0,"early"=60,"mid"=140,"late"=280)) {
  pop_df = get_muller_pop(x, mutations=mutations, label=label, map_tp_time=timepoints_to_int)
  highlight = get_highlight(x, min_frac, highlight)
  if (mutations) highlight = x %>%
    get_unique_muts_labels(clusters=highlight, label=label)
  color_palette = highlight_palette(x, highlight, label)

  p = pop_df %>% dplyr::filter(Identity %in% highlight) %>%
    ggplot(aes(x=Generation, y=Population, color=Identity)) +
    geom_point(alpha=.3) + my_ggplot_theme() + ylab("")

  for (cl in highlight) { p = exp_fit_util(p, pop_df, cl) }

  if (facet)
    p = p + scale_color_manual(values=color_palette[highlight]) +
    facet_wrap(Identity~Lineage, nrow=x$K) + xlab("Time") + labs(color="Clusters") else
      p = p + scale_color_manual(values=color_palette[highlight]) +
    facet_wrap(~Lineage) + xlab("Time") + labs(color="Clusters")
  return(p)
}


exp_fit_util = function(p, pop_df, cl) {
  lineages = pop_df %>% dplyr::filter(Identity == cl) %>% dplyr::pull(Lineage) %>% unique()
  exp_df = data.frame()
  for(ll in lineages) {
    lm_a = (pop_df %>% dplyr::filter(Identity==cl, Lineage==ll))$lm_a %>% unique()
    lm_r = (pop_df %>% dplyr::filter(Identity==cl, Lineage==ll))$lm_r %>% unique()

    xx = 1:max(pop_df$Generation)
    yy = exp(lm_a)*exp(lm_r*xx)
    exp_df = rbind(exp_df, data.frame(x=xx, y=yy, Identity=cl, Lineage=ll))
  }

  return(p + geom_line(data=exp_df, aes(x=x, y=y, color=Identity)))
}


#' Exponential rate
#'
#' @description Function to visualize the growth coefficients for each clone and lineage.
#'
#' @param x a mvnmm object.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param facet Boolean. If set to \code{TRUE}, the plot will be faceted against the clusters.
#' @param mutations Boolean. If set to \code{TRUE}, the growth will be visualize for each cluster of mutations.
#' @param label a character corresponding to the label of the run to visualize.
#' @param timepoints_to_int a list to map each \code{timepoint} value to an integer.
#'
#' @examples
#' if (FALSE) plot_exp_rate(x)
#'
#' @import ggplot2
#'
#' @export plot_exp_rate

plot_exp_rate = function(x, highlight=c(), min_frac=0, mutations=F, label="",
                         timepoints_to_int=list("init"=0,"early"=60,"mid"=140,"late"=280)) {

  highlight = get_highlight(x, min_frac, highlight, mutations=mutations, label=label) %>%
    color_palette = highlight_palette(x, highlight, label)

    highlight = highlight %>% stringr::str_replace("C_||C", "")
    names(color_palette) = names(color_palette) %>% stringr::str_replace("C_||C", "")

    p = x %>%
      get_muller_pop(map_tp_time=timepoints_to_int, mutations=mutations, label=label) %>%
      dplyr::mutate(Identity=Identity %>% stringr::str_replace("C_||C","")) %>%
      dplyr::filter(Identity %in% highlight) %>%
      dplyr::arrange(lm_r) %>%
      ggplot(aes(x=Identity, y=lm_r, ymax=lm_r, ymin=0, color=Identity)) +
      geom_linerange(alpha=.3) +
      geom_point(alpha=.3) +
      facet_wrap(~Lineage) +
      scale_color_manual(values=color_palette[highlight]) +
      my_ggplot_theme() +
      theme(aspect.ratio=1) +
      xlab("Clusters") + ylab("Exp rate") + labs(color="Clusters")
    return(p)
}
