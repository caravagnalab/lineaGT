plot_phylogeny = function(x, min_frac=0, highlight=c(), vaf.df=NULL, label="") {

  if (purrr::is_empty(highlight)) highlight = x %>% get_unique_labels()
  clusters_joined = intersect(select_relevant_clusters(x, min_frac), highlight)

  if (label!="") nn = paste("trees", label, sep=".")

  if (!"trees" %in% names(x)) {
    x = fit_phylogenies(x, vaf.df=vaf.df, min_frac=min_frac, highlight=highlight, label=label)
  }

  tree_plots = list()
  for (cluster in clusters_joined) {
    if (!purrr::is_empty(x$trees[[cluster]])) {
      tree = x[[nn]][[cluster]][[1]]
      tree_plots[[cluster]] = ctree:::plot.ctree(tree) + labs(title=cluster)
    }
  }

  return(tree_plots)
}



