plot_phylogeny = function(x, min_frac=0, highlight=c(), label="") {
  clusters_joined = get_highlight(x, min_frac, highlight)

  trees = get_trees(x, label)
  tree_plots = list()
  for (cluster in clusters_joined) {
    if (!purrr::is_empty(trees[[cluster]])) {
      tree = trees[[cluster]] %>% get_best_scores()

      for (tt in 1:length(tree))
        tree_plots[[cluster]][[tt]] = ctree:::plot.ctree(tree[[tt]]) + labs(title=cluster)
    }
  }

  return(tree_plots[[cluster]] %>% patchwork::wrap_plots())
}


get_best_scores = function(trees) {
  ret_trees = list()
  ret_trees[[1]] = trees[[1]]
  best_score = trees[[1]]$score
  for (idx in 2:length(trees)) {
    if (abs(trees[[idx]]$score - best_score) < 0.05)
      ret_trees[[idx]] = trees[[idx]]
  }
  return(ret_trees)
}

