plot_phylogeny = function(x, min_frac=0, highlight=c(), vaf.df=NULL) {

  if (purrr::is_empty(highlight)) highlight = x %>% get_unique_labels()
  clusters_joined = intersect(select_relevant_clusters(x, min_frac), highlight)

  if (!"trees" %in% names(x)) {
    x = fit_phylogenies(x, vaf.df=vaf.df, min_frac=min_frac, highlight=highlight)
  }

  tree_plots = list()
  for (cluster in clusters_joined) {
    if (!purrr::is_empty(x$trees[[cluster]])) {
      tree = x$trees[[cluster]][[1]]
      tree_plots[[cluster]] = ctree:::plot.ctree(tree) + labs(title=cluster)
    }
  }

  return(tree_plots)
}



