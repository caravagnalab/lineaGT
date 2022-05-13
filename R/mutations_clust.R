run_viber = function(x, vaf.df, min_frac=0, highlight=list()) {
  x = add_vaf(x, vaf.df)
  viber_input = get_input_viber(vaf.df, x)

  if (purrr::is_empty(highlight)) highlight = viber_input$vaf.df$labels %>% unique() %>% droplevels()
  clusters_joined = intersect(select_relevant_clusters(x, min_frac), highlight)

  joined = data.frame()
  fit_all = list()
  for (cluster in clusters_joined) {
    fit_k = fit_cluster_viber(viber_input, cluster=cluster)
    joined = rbind(joined, fit_k$df)
    fit_all[[cluster]] = fit_k$fit
  }

  x$viber_run = fit_all
  theta = get_binomial_theta(x)
  vaf.df = joined %>%
    mutate(labels_mut=paste(labels, labels_viber, sep=".")) %>%
    wide_to_long_muts()
  x$vaf.dataframe = dplyr::inner_join(vaf.df, theta, by=c("labels_mut","labels", "timepoints","lineage"))
  x$color_palette = c(x$color_palette, get_colors(list_lab=get_unique_muts_labels(x)))
  return(x)
}


fit_cluster_viber = function(viber_input, cluster) {
  viber_df_k = list("successes"=viber_input$successes %>% filter(labels==cluster) %>% dplyr::select(-labels),
                    "trials"=viber_input$trials %>% filter(labels==cluster) %>% dplyr::select(-labels),
                    "vaf.df"=viber_input$vaf.df %>% filter(labels==cluster))
  k = viber_df_k$successes %>% nrow
  fit_viber = list()
  try(expr = {
    fit_viber = VIBER::variational_fit(viber_df_k$successes, viber_df_k$trials, K=k)
    if (fit_viber$K * 0.01 < 1) pi_cutoff = 0.005 else pi_cutoff = 0.01
    fit_viber = VIBER::choose_clusters(fit_viber, binomial_cutoff=0,
                                       dimensions_cutoff=0, pi_cutoff=pi_cutoff, re_assign=T)

    labels = fit_viber$labels$cluster.Binomial
    viber_df_k$vaf.df$labels_viber = labels
    viber_df_k$vaf.df$pi_viber = fit_viber$pi_k[labels] %>% as.vector() }, silent = T)

  try(expr = {
    if (purrr::is_empty(fit_viber)) {
      viber_df_k$vaf.df$labels_viber = "N"
      viber_df_k$vaf.df$pi_viber = 0
    }
  }, silent = T )

  return(list("df"=viber_df_k$vaf.df, "fit"=fit_viber))
}


