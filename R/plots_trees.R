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
#' @export plot_phylogeny

plot_phylogeny = function(x, score_diff=1, show_best=1, min_frac=0, highlight=c(), label="") {
  clusters_joined = get_highlight(x, min_frac, highlight)

  trees = get_trees(x, label)
  tree_plots = list()
  for (cluster in clusters_joined) {
    if (!purrr::is_empty(trees[[cluster]])) {
      tree = trees[[cluster]] %>% get_best_scores(show_best=show_best, score_diff=score_diff)

      for (tt in 1:length(tree))
        tree_plots[[cluster]][[tt]] = ctree:::plot.ctree(tree[[tt]]) + labs(title=cluster)
    }
  }

  return(tree_plots[[cluster]] %>% patchwork::wrap_plots())
}


get_best_scores = function(trees, show_best=0, score_diff=1) {
  if (show_best == 1 && score_diff == 1)
    return(trees[1])

  if (show_best == 0 || score_diff < 1)
    show_best = length(trees)

  ret_trees = list()
  ret_trees[[1]] = trees[[1]]
  best_score = trees[[1]]$score

  max_n = min(show_best, length(trees))
  for (idx in 2:max_n) {
    if (abs(trees[[idx]]$score - best_score) < score_diff)
      ret_trees[[idx]] = trees[[idx]]
    # ret_trees[[idx]] = trees[[idx]]
  }

  return(ret_trees)
}

