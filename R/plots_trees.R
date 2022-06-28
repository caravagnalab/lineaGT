#' Clonal evolution trees
#'
#' @param x add
#' @param score_diff add
#' @param show_best add
#' @param min_frac add
#' @param highlight add
#' @param label add
#'
#' @return
#'
#' @importFrom ggraph create_layout ggraph geom_edge_link geom_node_point geom_node_text circle
#' @importFrom RColorBrewer brewer.pal
#'
#' @export plot_phylogeny

plot_phylogeny = function(x, score_diff=1, show_best=1, min_frac=0, highlight=c(), label="") {
  clusters_joined = get_highlight(x, min_frac, highlight)

  trees = get_trees(x, label)
  tree_plots = list()
  for (cluster in clusters_joined) {
    if (!purrr::is_empty(trees[[cluster]])) {
      tree = trees[[cluster]] %>% get_best_scores(show_best=show_best, score_diff=score_diff)

      color_palette = x %>% get_color_palette()
      color_palette = color_palette[x %>% get_unique_muts_labels(clusters=cluster)]
      names(color_palette) = names(color_palette) %>%
        str_replace_all(cluster, "") %>%
        str_replace_all("[.]", "")

      cluster_plots = list()
      for (tt in 1:length(tree))
        cluster_plots[[tt]] = plot_ctree_mod(tree[[tt]],
                                             node_palette=color_palette,
                                             cluster_id=cluster)
      tree_plots[[cluster]] = patchwork::wrap_plots(cluster_plots, guides="collect") & theme(legend.position = "bottom")
    }
  }

  return(tree_plots)
}



plot_ctree_mod = function (x.tree,
                           cluster_id="",
                           node_palette=c(),
                           tree_layout="tree", ...) {
  if (purrr::is_empty(node_palette))
    node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))

  tree = x.tree
  cex = 1
  tb_tree = tree$tb_adj_mat
  clusters = tree$CCF %>% dplyr::pull(cluster) %>% unique()
  node_palette[setdiff(clusters, names(node_palette))] = "gainsboro"

  layout = ggraph::create_layout(tb_tree, layout=tree_layout)
  return(
    ggraph::ggraph(layout) +
      ggraph::geom_edge_link(arrow=arrow(length=unit(2 * cex, "mm")),
                             end_cap=ggraph::circle(5 * cex, "mm"),
                             start_cap=ggraph::circle(5 * cex, "mm")) +
      ggraph::geom_node_point(aes(colour=cluster, size=nMuts), alpha=.6, na.rm=TRUE) +
      ggraph::geom_node_text(aes(label=cluster), colour="black", vjust=0.4) +
      coord_cartesian(clip="off") +
      scale_color_manual(values=node_palette) +
      scale_size(range=c(3, 10) * cex) +
      guides(color=FALSE, size=guide_legend("Clone size", nrow=1)) +
      labs(title=paste(cluster_id),
           subtitle=paste0("Scores ", format(tree$score, scientific=T), ".")) +
      theme_void(base_size=10*cex) +
      theme(legend.position="bottom", text=element_text(size=9))
    )
}


get_best_scores = function(trees, show_best=0, score_diff=1) {
  if ((show_best == 1 && score_diff == 1) || length(trees)==1)
    return(trees[1])

  if (show_best == 0 || score_diff < 1)
    show_best = min(length(trees), 3)

  ret_trees = list()
  ret_trees[[1]] = trees[[1]]
  best_score = trees[[1]]$score

  max_n = min(show_best, length(trees))
  for (idx in 2:max_n) {
    if (abs(trees[[idx]]$score - best_score) < score_diff)
      ret_trees[[idx]] = trees[[idx]]
  }

  return(ret_trees)
}

