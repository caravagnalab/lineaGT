#' VAF 2D scatterplot
#'
#' @description Function to plot the VAFs of the mutations one timepoint against the other
#'
#' @param x a \code{mvnmm} object.
#' @param min_frac value in \code{[0,1]} to show only the clusters with a minimum CCF of
#' \code{min_frac} in at least one timepoint.
#' @param highlight a list of labels ID to show. All the clusters in the list are shown,
#' even if the CCF is lower than \code{min_frac}.
#' @param wrap a Boolean to wrap the scatterplots of multiple lineages in a unique plot.
#' @param label a character corresponding to the label of the fit to visualize, if more than one.
#'
#' @return list of VAF scatterplots.
#'
#' @examples
#' if (FALSE) plot_vaf(x, min_frac=0.1)
#'
#' @import ggplot2
#' @importFrom dplyr select starts_with
#' @importFrom purrr is_empty
#' @importFrom patchwork wrap_plots
#'
#' @export plot_vaf

plot_vaf = function(x, min_frac=0, highlight=c(), wrap=T, label="") {
  if (label != "")
    dataframe = x %>% add_lineage_vaf(label=label)
  else
    dataframe = x %>% get_vaf_dataframe()

  dataframe = dataframe %>%
    dplyr::mutate(theta=theta*100) %>%
    dplyr::select(-contains("ref"), -contains("dp"), -contains("alt")) %>%
    tidyr::pivot_wider(names_from=c("timepoints"),
                       names_sep=".",
                       values_from=c("vaf","theta"))

  highlight = get_highlight(x, min_frac, highlight, mutations=T, label=label)
  highlight_v = get_unique_muts_labels(x, highlight, label=label)
  color_palette = highlight_palette(x, highlight, label=label)

  combinations = get_pairs(dataframe, columns=dataframe %>%
                             dplyr::select(dplyr::starts_with("vaf")) %>%
                             colnames)

  p = list()
  for (t1_t2 in combinations$pair_name) {
    xy.vaf = strsplit(t1_t2, ":")[[1]]
    xy.theta = str_replace_all(xy.vaf, "vaf.", "theta.")

    df = dataframe %>%
      filter(labels %in% highlight)

    if (nrow(df) > 0) { p[[t1_t2]] = plot_vaf_2D(df, xy.vaf, xy.theta, color_palette[highlight_v]) }
  }

  if (wrap) return(p %>% patchwork::wrap_plots(guides="collect"))
  return(p)
}



plot_vaf_2D = function(dataframe, dims.vaf, dims.theta, color_palette) {
  if (nrow(dataframe) == 1)
    return(message("Only one mutation"))
  pl = dataframe %>%
    dplyr::select(starts_with("vaf"), starts_with("theta"), labels, lineage, labels_mut, mutation, pi_viber) %>%
    ggplot() +
    geom_point(aes_string(x=dims.vaf[1], y=dims.vaf[2], color="labels_mut"), alpha=.5, size=.7) +
    geom_point(aes_string(x=dims.theta[1], y=dims.theta[2], color="labels_mut"), shape=15, inherit.aes=F, size=1.5) +
    facet_grid(lineage ~ labels) +
    scale_color_manual(values=color_palette) +
    xlab(split_to_camelcase(dims.vaf[1])) +
    ylab(split_to_camelcase(dims.vaf[2])) +
    ylim(0, 100) +
    xlim(0, 100) +
    labs(color = "Clusters") +
    my_ggplot_theme()

  return(pl)
}



#' VAF over time
#'
#' @description Function to plot the VAFs of the mutations along time
#'
#' @param x a \code{mvnmm} object.
#' @param min_frac value in \code{[0,1]} to show only the clusters with a minimum CCF of
#' \code{min_frac} in at least one timepoint.
#' @param highlight a list of labels ID to show. All the clusters in the list are shown,
#' even if the CCF is lower than \code{min_frac}.
#' @param timepoints_to_int a list to map each \code{timepoint} value to an integer.
#' @param label a character corresponding to the label of the fit to visualize, if more than one.
#'
#' @return a \code{ggplot} object.
#'
#' @examples
#' if (FALSE) plot_vaf_time(x, min_frac=0.1)
#'
#' @import ggplot2
#' @importFrom dplyr select starts_with
#' @importFrom purrr is_empty
#' @importFrom patchwork wrap_plots
#'
#' @export plot_vaf_time

plot_vaf_time = function(x, min_frac=0, highlight=c(),
                         timepoints_to_int=list("init"=0,"early"=60,"mid"=140,"late"=280),
                         label="") {

  highlight.c = get_highlight(x, min_frac, highlight, mutations=F)
  highlight.m = get_unique_muts_labels(x, clusters=highlight.c, label=label)

  vaf.df = x %>% add_lineage_vaf(label=label) %>%
    filter(labels %in% highlight.c) %>%
    mutate(timepoints=timepoints_to_int[timepoints]) %>%
    mutate(timepoints=unlist(timepoints)) %>%
    mutate(timepoints=as.numeric(timepoints))

  color_palette = get_color_palette(x,label=label)[highlight.m]
  color_palette[paste(vaf.df %>% dplyr::pull(timepoints) %>% unique())] = "black"

  return(
    vaf.df %>%
      ggplot() +
      geom_point(aes(x=timepoints, y=vaf, color=labels_mut), alpha=.5, size=.7) +
      geom_line(aes(x=timepoints, y=vaf, color=labels_mut, group=mutation),
                alpha=.6, size=.5, linetype="solid") +
      geom_point(aes(x=timepoints, y=theta*100, color=labels_mut), shape=15, size=1.5) +
      geom_vline(xintercept=vaf.df$timepoints %>% unique(), linetype="dashed", size=.4, alpha=.5) +
      facet_wrap(labels~lineage, ncol=x %>% get_lineages() %>% length()) +
      scale_color_manual(values=color_palette, breaks=highlight.m) +
      ylab("VAF") + xlab("Time") + labs(color="Clusters") +
      ylim(0,100) + my_ggplot_theme()
  )
}



add_lineage_vaf = function(x, label="") {
  if (label == "")
    return(x %>% get_vaf_dataframe())

  ll_used = x %>% get_vaf_dataframe(label=label) %>% dplyr::pull(lineage) %>% unique()
  ll_notused = (x %>% get_lineages())[! (x %>% get_lineages()) %in% ll_used]

  vaf.1 = x %>%
    get_vaf_dataframe() %>%
    filter(lineage %in% ll_notused) %>%
    dplyr::select(-labels_viber, -pi_viber, -labels_mut, -theta) %>%
    long_to_wide_muts()

  vaf.2 = x %>%
    get_vaf_dataframe(label=label) %>%
    long_to_wide_muts()

  return(
    inner_join(vaf.1, vaf.2, by=c("IS", "mutation", "labels", "labels_init")) %>%
      wide_to_long_muts()
  )
}
