# `label` is to give a label to the run
# `lineages` is to select on which lineages to run the inference,
# otherwise it will be performed on all of them
run_viber = function(x, vaf.df=NULL, min_frac=0, highlight=list(), do_filter=FALSE,
                     infer_phylo=TRUE, lineages=c(), label="") {

  if (is.null(vaf.df) && !"vaf.dataframe" %in% names(x))
    message("A dataframe with the mutations is required!")
  else if (is.null(vaf.df)) vaf.df = x %>% get_vaf_dataframe(label)

  clusters_joined = get_highlight(x, min_frac, highlight)
  x = x %>%
    annotate_vaf_df(vaf.df=vaf.df, min_frac=min_frac, label=label) %>%  # add cluster to each mutation
    check_dp(label=label)  # check for too low values

  viber_input = x %>% get_input_viber(lineages=lineages, label=label)  # get the input to run viber

  joined = data.frame(); fit_all = list()
  if (infer_phylo) trees = list()

  for (cluster in clusters_joined) {
    fit_k = fit_cluster_viber(viber_input,
                              cluster=cluster,
                              infer_phylo=infer_phylo)

    joined = rbind(joined, fit_k$df)
    fit_all[[cluster]] = fit_k$fit
    if (infer_phylo)
      trees[[cluster]] = fit_k$tree
  }

  x = add_viber_run(x, viber_run=fit_all, label=label)

  vaf.df = add_theta_to_vaf(x, joined, label=label)
  x = add_vaf(x, vaf.df=vaf.df, label=label)

  color_palette = update_color_palette(x, clusters_joined, label=label)
  x = add_color_palette(x, color_palette, label=label)

  x = add_phylo(x, trees, label=label)

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

