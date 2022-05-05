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
#' @import ggplot2
#' @importFrom dplyr select starts_with
#' @importFrom purrr is_empty
#'
#' @export plot_vaf


plot_vaf = function(x, min_frac=0, highlight=c()) {
  dataframe = reshape_vaf_dataframe_long(x)

  if (purrr::is_empty(highlight)) highlight = select_relevant_clusters(x, min_frac)
  highlight_v = get_viber_clusters(x, highlight)
  color_palette = highlight_palette(x$color_palette, c(highlight, highlight_v))

  combinations = get_pairs(dataframe, columns=dataframe %>% dplyr::select(dplyr::starts_with("vaf")) %>% colnames)
  theta = reshape_theta_long(x)

  p = list()
  for (t1_t2 in combinations$pair_name) {
    xy = strsplit(t1_t2, ":")[[1]]

    df = dataframe %>% filter(labels %in% highlight)
    tt = theta %>% filter(labels %in% highlight)
    if (nrow(df) > 0) { p[[t1_t2]] = plot_vaf_2D(df, tt, xy[1], xy[2], color_palette[highlight_v]) }
  }
  return(p)
}



plot_vaf_2D = function(dataframe,
                       theta,
                       dim1,
                       dim2,
                       color_palette) {
  pl = dataframe %>%
    dplyr::select(starts_with("vaf"),
                  labels,
                  lineage,
                  labels_mut,
                  mutation,
                  pi_viber) %>%
    ggplot() +
    geom_point(aes_string(x = dim1, y = dim2, color = "labels_mut"), alpha =
                            .5) +
    geom_point(
      data = theta,
      aes_string(x = dim1, y = dim2, color = "labels_mut"),
      shape = 15,
      inherit.aes = F
    ) +
    facet_grid(lineage ~ labels) +
    scale_color_manual(values = color_palette) +
    my_ggplot_theme() +
    xlab(split_to_camelcase(dim1)) +
    ylab(split_to_camelcase(dim2)) +
    ylim(0, 100) + xlim(0, 100) + labs(color = "Clusters")

  return(pl)
}

