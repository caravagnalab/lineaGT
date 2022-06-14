plot_phylogeny = function(x, min_frac=0, highlight=c(), label="") {
  clusters_joined = get_highlight(x, min_frac, highlight)

  trees = get_trees(x, label)
  tree_plots = list()
  for (cluster in clusters_joined) {
    if (!purrr::is_empty(trees[[cluster]])) {
      tree = trees[[cluster]][[1]]
      tree_plots[[cluster]] = ctree:::plot.ctree(tree) + labs(title=cluster)
    }
  }

  return(tree_plots)
}

