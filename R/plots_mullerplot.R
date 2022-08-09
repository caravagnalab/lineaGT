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
#' @param tree_score add
#' @param legend.pos add
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

plot_mullerplot = function(x,
                           which="frac",
                           highlight=c(),
                           min_frac=0,
                           timepoints_to_int=c(),
                           mutations=F,
                           single_clone=T,
                           tree_score=1,
                           legend.pos="right",
                           wrap=T) {

  timepoints_to_int = map_timepoints_int(x, timepoints_to_int)

  highlight.cov = get_highlight(x, min_frac=min_frac, highlight=highlight)
  highlight = get_highlight(x, min_frac, highlight.cov, mutations=mutations)
  color_palette = highlight_palette(x, highlight)

  lvls = c("P", get_unique_muts_labels(x), get_unique_labels(x))

  pop_df = get_muller_pop(x, mutations=mutations, timepoints_to_int=timepoints_to_int) %>%
    dplyr::select(-Population, -Frequency, -Parent, -theta_binom, -dplyr::contains("Pop.subcl")) %>%
    dplyr::rename(Population=Pop.plot) %>%

    dplyr::mutate(Identity=factor(Identity, levels=lvls)) %>%
    dplyr::arrange(Identity, Generation, Lineage)

  edges_df = get_muller_edges(x,
                              mutations=mutations,
                              tree_score=tree_score) %>%
    dplyr::mutate(Parent=factor(Parent, levels=lvls)) %>%
    dplyr::arrange(Parent) #%>% dplyr::mutate(Parent=as.character(Parent))


  if (single_clone && mutations) {
    pop_df = pop_df %>% filter_muller_df(highlight=highlight.cov)
    edges_df = edges_df %>% filter_muller_df(highlight=highlight.cov)
  }

  timepoints = x %>% get_dimensions()
  lineages = x %>% get_lineages()

  mullerdf = data.frame()
  for (ll in lineages) {
    tp = timepoints[grep(pattern=ll, x=timepoints)]
    if (length(tp) != 0) {
      pop_ll = pop_df %>% dplyr::filter(Lineage==ll)
      mullerdf_ll = ggmuller::get_Muller_df(edges_df, pop_ll) %>%
        dplyr::mutate(Lineage=ll)

      mullerdf = rbind(mullerdf, mullerdf_ll)
    }
  }

  return(
    mullerplot_util(mullerdf %>% filter(Generation %in% (timepoints_to_int %>% unlist())),
                    which=which,
                    highlight=highlight,
                    color_palette=color_palette,
                    legend.pos=legend.pos)
  )
}


mullerplot_util = function(mullerdf, which, color_palette, highlight, legend.pos="right") {
  if (which == "fitness") {
    exp_limits = c(min(mullerdf$lm_r), max(mullerdf$lm_r))

    return(
      mullerdf %>% ggplot() +
        geom_area(aes_string(x="Generation", y="Frequency", group="Group_id", fill="lm_r")) +
        geom_vline(xintercept=mullerdf$Generation %>% unique(), linetype="dashed") +
        guides(linetype="none", color="none") +
        facet_wrap(~Lineage, nrow=1) +
        xlab("Time") +
        labs(fill="Exp rate") +
        my_ggplot_theme(legend.pos=legend.pos) +
        scale_fill_gradient2(mid="white", low="blue", high="red", limits=exp_limits, na.value="#FFFFFF00")
    )
  }

  if (which == "frac")
    y = "Frequency"
  else if (which == "pop") {
    y = "Population"
    mullerdf = mullerdf %>% pop_df_add_empty()
  }

  return(
    mullerdf %>%
      dplyr::mutate(Identity=factor(Identity, levels=lvls)) %>%
      dplyr::arrange(Identity, Group_id) %>%
      ggplot() +
      geom_area(aes_string(x="Generation", y=y, group="Group_id", fill="Identity", colour="Identity")) +
      geom_vline(xintercept=mullerdf$Generation %>% unique(), linetype="dashed") +
      guides(linetype="none", color="none") +
      facet_wrap(~Lineage, nrow=1) +
      scale_fill_manual(name="Clusters", values=color_palette, na.value="#FFFFFF00", breaks=highlight) +
      scale_color_manual(values=color_palette, na.value="#FFFFFF00", breaks=highlight) +
      xlab("Time") +
      my_ggplot_theme(legend.pos=legend.pos)
  )

}

