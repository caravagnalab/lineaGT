#' VAF scatterplot
#'
#' @description Function to plot the VAFs of the mutations one timepoint against the other
#'
#' @param x a \code{mvnmm} object.
#' @param min_fracc value in \code{[0,1]} to show only the clusters with a minimum CCF of
#' \code{min_frac} in at least one timepoint.
#' @param highlight a list of labels ID to show. All the clusters in the list are shown,
#' even if the CCF is lower than \code{min_frac}.
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


plot_vaf = function(x, min_frac=0, highlight=c(), label="", wrap=T) {
  if (label != "")
    dataframe = x %>% add_lineage_vaf(label=label)
  else
    dataframe = x %>% get_vaf_dataframe()

  dataframe = dataframe %>%
    dplyr::select(-contains("ref"), -contains("dp"), -contains("alt"), -contains("theta")) %>%
    tidyr::pivot_wider(names_from=c("timepoints"), names_sep=".", values_from=c("vaf"), names_prefix="vaf.")

  highlight = get_highlight(x, min_frac, highlight, mutations=T, label=label)
  highlight_v = get_unique_muts_labels(x, highlight, label=label)
  color_palette = highlight_palette(x, highlight, label=label)

  combinations = get_pairs(dataframe, columns=dataframe %>%
                             dplyr::select(dplyr::starts_with("vaf")) %>%
                             colnames)
  theta = x %>% get_binomial_theta(label=label) %>%
    tidyr::pivot_wider(names_from=c("timepoints"),
                       values_from=c("theta"),
                       names_prefix="vaf.",
                       names_sep=".")

  p = list()
  for (t1_t2 in combinations$pair_name) {
    xy = strsplit(t1_t2, ":")[[1]]

    df = dataframe %>% filter(labels %in% highlight)
    tt = theta %>% filter(labels %in% highlight)
    if (nrow(df) > 0) { p[[t1_t2]] = plot_vaf_2D(df, tt, xy[1], xy[2], color_palette[highlight_v]) }
  }

  if (wrap) return(p %>% patchwork::wrap_plots(guides="collect"))
  return(p)
}



plot_vaf_2D = function(dataframe, theta, dim1, dim2, color_palette) {
  pl = dataframe %>%
    dplyr::select(starts_with("vaf"), labels, lineage, labels_mut, mutation, pi_viber) %>%
    ggplot() +
    geom_point(aes_string(x=dim1, y=dim2, color="labels_mut"), alpha=.5, size=.7) +
    geom_point(data=theta, aes_string(x=dim1, y=dim2, color="labels_mut"), shape=15, inherit.aes=F, size=1.5) +
    facet_grid(lineage ~ labels) +
    scale_color_manual(values=color_palette) +
    xlab(split_to_camelcase(dim1)) +
    ylab(split_to_camelcase(dim2)) +
    ylim(0, 100) +
    xlim(0, 100) +
    labs(color = "Clusters") +
    my_ggplot_theme()

  return(pl)
}



add_lineage_vaf = function(x, label="") {
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

