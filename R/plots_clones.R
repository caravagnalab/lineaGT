#' Muller plot
#'
#' @description
#'
#' @import ggplot2
#' @import ggmuller
#' @importFrom patchwork wrap_plots
#'
#' @export


plot_mullerplot = function(obj, which="ccf", highlight=c(), min_ccf=0, legend.pos="right", wrap=F, viber=F) {
  if (viber) {
    theta = get_binomial_theta(obj)
    pop_df = get_muller_pop(obj, means=theta)
    edges_df = get_muller_edges(obj, labels=get_unique_viber_labels(obj))

    obj$vaf_dataframe = get_vaf_dataframe(obj) %>% mutate(labels_mut=paste(labels,labels_viber,sep="."))
    if (!purrr::is_empty(highlight))
      highlight = get_unique_viber_labels(obj)[grepl(highlight, get_unique_viber_labels(obj))] else
        highlight = select_relevant_clusters(obj, min_ccf, theta)
    color_palette = highlight_palette(obj$color_palette, highlight)

  } else {
    if (purrr::is_empty(highlight)) highlight = select_relevant_clusters(obj, min_ccf)
    color_palette = highlight_palette(obj$color_palette, highlight)
    pop_df = get_muller_pop(obj)
    edges_df = get_muller_edges(obj)
  }

  timepoints = obj$dimensions
  lineages = obj$lineages
  exp_limits = c(min(pop_df$lm_r), max(pop_df$lm_r))

  plot_list = list()
  for (ll in lineages) {
    tp = timepoints[grep(pattern=ll, x=timepoints)]
    if (length(tp) != 0) {
      pop_ll = pop_df %>% filter(Lineage==ll)
      mullerdf_ll = ggmuller::get_Muller_df(edges_df, pop_ll)

      if (which == "ccf" || which == "")
        plot_list[[ll]] = mullerplot_util(mullerdf_ll, y="Frequency", fill="Identity", lineage=ll,
                                          color_palette=color_palette[highlight], legend.pos=legend.pos)
      if (which == "pop" || which == "")
        plot_list[[ll]] = mullerplot_util(mullerdf_ll %>% ggmuller::add_empty_pop(), y="Population",
                                          fill="Identity", color_palette=color_palette[highlight],
                                          lineage=ll, legend.pos=legend.pos)
      if (which == "fitness")
        plot_list[[ll]] = mullerplot_util(mullerdf_ll, y="Frequency", fill="lm_r",
                                          color_palette=color_palette[highlight], lineage=ll,
                                          legend.pos=legend.pos, exp_limits=exp_limits)
    }
  }

  if (wrap) return(patchwork::wrap_plots(plot_list, guides="collect"))
  return(plot_list)
}


mullerplot_util = function(mullerdf, y, fill, lineage, color_palette, legend.pos="right", exp_limits=NULL) {
  if (fill=="Identity")
    pl = mullerdf %>% ggplot() +
      geom_area(aes_string(x="Generation", y=y, group="Group_id", fill="Identity", colour="Identity"), alpha=.9) +
      guides(linetype=FALSE, color=FALSE) +
      xlab("Time") + labs(title=split_to_camelcase(lineage)) +
      my_ggplot_theme(legend.pos=legend.pos) +
      scale_fill_manual(name="Clusters", values=color_palette, na.value="transparent") +
      scale_color_manual(values=color_palette, na.value="transparent")

  if (fill == "lm_r")
    pl = mullerdf %>% ggplot() +
      geom_area(aes_string(x="Generation", y=y, group="Group_id", fill="lm_r"), alpha=.9) +
      guides(linetype=FALSE, color=FALSE) +
      xlab("Time") + labs(title=split_to_camelcase(lineage), fill="Exp rate") +
      my_ggplot_theme(legend.pos=legend.pos) +
      scale_fill_gradient2(mid="white", low="blue", high="red", limits=exp_limits)

  return(pl)
}


#' Exponential fitting
#'
#' @description
#'
#' @import ggplot2
#'
#' @export


plot_exp_fit = function(obj, highlight=c(), min_ccf=0, facet=F, viber=F) {
  if (viber) {
    theta = get_binomial_theta(obj)
    pop_df = get_muller_pop(obj, means=theta)

    obj$vaf_dataframe = get_vaf_dataframe(obj) %>% mutate(labels_mut=paste(labels,labels_viber,sep="."))
    if (purrr::is_empty(highlight)) highlight_c = select_relevant_clusters(obj, min_ccf) else highlight_c = highlight
    highlight_v = (dataframe %>% filter(labels %in% highlight_c))$labels_mut %>% unique()
    color_palette = highlight_palette(obj$color_palette, c(highlight_c,highlight_v))

  } else {
    if (purrr::is_empty(highlight)) highlight = select_relevant_clusters(obj, min_ccf)
    color_palette = highlight_palette(obj$color_palette, highlight)
    pop_df = get_muller_pop(obj)
  }

  p = pop_df %>% filter(Identity %in% highlight) %>%
    ggplot(aes(x=Generation, y=Population, color=Identity)) +
    geom_point(alpha=.3) + my_ggplot_theme()

  for (cl in highlight_c) { p = exp_fit_util(p, pop_df, cl) }

  if (facet)
    p = p + scale_color_manual(values=color_palette[highlight]) +
    facet_wrap(Identity~Lineage, nrow=obj$K) + xlab("Time") + labs(color="Clusters") else
      p = p + scale_color_manual(values=color_palette[highlight]) +
    facet_wrap(~Lineage) + xlab("Time") + labs(color="Clusters")
  return(p)
}


exp_fit_util = function(p, pop_df, cl) {
  exp_df = data.frame()
  for(ll in pop_df$Lineage %>% unique()) {
    lm_a = (pop_df %>% filter(Identity==cl, Lineage==ll))$lm_a %>% unique()
    lm_r = (pop_df %>% filter(Identity==cl, Lineage==ll))$lm_r %>% unique()

    xx = 1:max(pop_df$Generation)
    yy = exp(lm_a)*exp(lm_r*xx)
    exp_df = rbind(exp_df, data.frame(x=xx, y=yy, Identity=cl, Lineage=ll))
  }

  return(p + geom_line(data=exp_df, aes(x=x, y=y, color=Identity)))
}


#' Exponential rate
#'
#' @description
#'
#' @import ggplot2
#'
#' @export


plot_exp_rate = function(obj, highlight=c(), min_ccf=0) {
  if (purrr::is_empty(highlight)) highlight = select_relevant_clusters(obj, min_ccf)
  color_palette = highlight_palette(obj$color_palette, highlight)

  highlight = highlight %>% stringr::str_replace("C_", "")
  names(color_palette) = names(color_palette) %>% stringr::str_replace("C_", "")

  p = get_muller_pop(obj) %>% mutate(Identity=Identity %>% stringr::str_replace("C_","")) %>%
    filter(Identity%in%highlight) %>%
    dplyr::arrange(lm_r) %>%
    ggplot(aes(x=Identity, y=lm_r, ymax=lm_r, ymin=0, color=Identity)) +
    geom_linerange(alpha=.3) + geom_point(alpha=.3) + facet_wrap(~Lineage) +
    scale_color_manual(values=color_palette[highlight]) + my_ggplot_theme() +
    theme(aspect.ratio=1) + xlab("Clusters") + ylab("Exp rate") + labs(color="Clusters")
  return(p)
}
