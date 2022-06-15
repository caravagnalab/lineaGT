# function to infer the phylogenies on a fit object
fit_phylogenies = function(x, vaf.df=NULL, min_frac=0, highlight=list(), do_filter=FALSE,
                           label="", fit_viber=FALSE, lineages=c()) {

  if (is.null(vaf.df) && !"vaf.dataframe" %in% names(x))
    message("A dataframe with the mutations is required!")
  else if (is.null(vaf.df))
    vaf.df = x %>% get_vaf_dataframe(label)

  clusters_joined = get_highlight(x, min_frac, highlight)
  trees = list()
  if (!"viber_run" %in% names(x) || fit_viber)
    return(
      x %>% run_viber(vaf.df=vaf.df,
                      highlight=clusters_joined,
                      lineages=lineages,
                      label=label,
                      infer_phylo=TRUE)
    )

  viber_run_all = x %>% get_viber_run(label=label)

  for (cluster in clusters_joined) {
    viber_run = viber_run_all[[cluster]]
    tt = fit_trees(viber_run)
    trees[[cluster]] = tt
  }

  x = add_phylo(x, trees, label=label)

  return(x)
}


# to infer the tree on a single cluster
fit_trees = function(fit_viber) {
  if (length(fit_viber$labels$cluster.Binomial %>% unique) > 1)
    tree = VIBER::get_clone_trees(fit_viber)

  return(tree)
}

