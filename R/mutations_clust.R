run_viber = function(x, vaf.df, min_frac=0, highlight=list(), do_filter=FALSE, infer_phylo=TRUE) {
  clusters_joined = retrieve_clusters(x, min_frac, highlight)

  x = x %>%
    annotate_vaf_df(vaf.df=vaf.df, min_frac=min_frac) %>%  # add cluster to each mutation
    check_dp()  # check for too low values

  viber_input = x %>% get_input_viber()  # get the input to run viber

  joined = data.frame(); fit_all = list()
  if (infer_phylo) trees = list() # plots = list()

  for (cluster in clusters_joined) {
    fit_k = fit_cluster_viber(viber_input, cluster=cluster, infer_phylo=infer_phylo)
    joined = rbind(joined, fit_k$df)
    fit_all[[cluster]] = fit_k$fit

    if (infer_phylo) {
      trees[[cluster]] = fit_k$tree
    }
  }

  x$viber_run = fit_all
  x$trees = trees
  theta = get_binomial_theta(x)
  vaf.df = joined %>%
    dplyr::mutate(labels_mut=paste(labels, labels_viber, sep=".")) %>%
    wide_to_long_muts()
  x$vaf.dataframe = dplyr::inner_join(vaf.df, theta, by=c("labels_mut","labels", "timepoints","lineage"))
  x$color_palette = c(x$color_palette,
                      get_colors(x=x,
                                 list_lab=get_unique_muts_labels(x),
                                 color_palette=x$color_palette))

  return(x)
}


fit_cluster_viber = function(viber_input, cluster, infer_phylo=TRUE) {
  viber_df_k = list("successes"=viber_input$successes %>%
                      dplyr::filter(labels==cluster) %>%
                      dplyr::select(-labels),
                    "trials"=viber_input$trials %>%
                      dplyr::filter(labels==cluster) %>%
                      dplyr::select(-labels),
                    "vaf.df"=viber_input$vaf.df %>%
                      dplyr::filter(labels==cluster))
  k = viber_df_k$successes %>% nrow  # max n of clusters
  fit_viber = tree_joint = list()

  try(expr = {
    data_annotations = data.frame(gene=paste0("G", 1:k), driver=FALSE)
    data_annotations$driver[sample(1:nrow(data_annotations), 1)] = TRUE

    fit_viber = VIBER::variational_fit(viber_df_k$successes,
                                       viber_df_k$trials,
                                       K=k,
                                       data=data_annotations)
    if (fit_viber$K * 0.01 < 1) pi_cutoff = 0.005 else pi_cutoff = 0.01
    fit_viber = VIBER::choose_clusters(fit_viber,
                                       binomial_cutoff=0,
                                       dimensions_cutoff=0,
                                       pi_cutoff=pi_cutoff,
                                       re_assign=T)

    labels = fit_viber$labels$cluster.Binomial
    viber_df_k$vaf.df$labels_viber = labels
    viber_df_k$vaf.df$pi_viber = fit_viber$pi_k[labels] %>% as.vector()

    if (infer_phylo) {
      tt = fit_trees(fit_viber)
      tree_joint = tt
    }

  }, silent = T)

  try(expr = {
    if (purrr::is_empty(fit_viber)) {
      viber_df_k$vaf.df$labels_viber = "N"
      viber_df_k$vaf.df$pi_viber = 0
    }
  }, silent = T )

  return(list("df"=viber_df_k$vaf.df, "fit"=fit_viber, "tree"=tree_joint)) #, "plot"=plot_joint))
}


fit_phylogenies = function(x, vaf.df=NULL, min_frac=0, highlight=list(), do_filter=FALSE) {
  clusters_joined = retrieve_clusters(x, min_frac, highlight)

  trees = list()

  if (!"viber_run" %in% names(x)) {
    x = x %>% run_viber(vaf.df=vaf.df, highlight=clusters_joined)
  }

  for (cluster in clusters_joined) {
    viber_run = x$viber_run[[cluster]]
    tt = fit_trees(viber_run)
    trees[[cluster]] = tt
  }

  x$trees = trees

  return(x)
}


# to infer the tree on a single cluster
fit_trees = function(fit_viber) {
  tree = list() # plot = list()
  if (length(fit_viber$labels$cluster.Binomial %>% unique) > 1) {
    tree = VIBER::get_clone_trees(fit_viber)
  }

  return(tree) #, "plot"=plot))
}

