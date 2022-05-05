#' Functions to run VIBER for mutations clustering
#'
#' @description add
#'
#' @param x mvnmm object.
#' @param vaf_df a VAF dataframe, with required columns mutation,dp,ref,alt of depth, reference and alternative reads
#' @param min_frac numeric value in \code{[0,1]} representing the minimum abundance to consider a clone in the inference.
#' @return mvnmm object with added the VAF dataframe annotated with the found clusters and the VIBER fit for each cluster
#'
#' @export run_viber


run_viber = function(x, vaf_df, min_frac=0, highlight=list()) {
  x = add_vaf(x, vaf_df)

  try(expr = { vaf_df = vaf_df %>% filter(mutation %in% x$keep_mut) }, silent = T)
  viber_input = get_input_viber(vaf_df)

  if (purrr::is_empty(highlight)) highlight = viber_input$vaf_df$labels %>% unique() %>% droplevels()
  clusters_joined = intersect(select_relevant_clusters(x, min_frac), highlight)

  joined = data.frame()
  fit_all = list()
  for (cluster in clusters_joined) {
    fit_k = fit_cluster_viber(viber_input, cluster=cluster)
    joined = rbind(joined, fit_k$df)
    fit_all[[cluster]] = fit_k$fit
  }
  x$viber_run = fit_all
  x$vaf_all = x$vaf_dataframe
  x$vaf_dataframe = joined %>% mutate(labels_mut=paste(labels, labels_viber, sep="."))
  x$color_palette = c(x$color_palette, get_colors(list_lab=get_unique_viber_labels(x)))
  return(x)
}


fit_cluster_viber = function(viber_input, cluster) {
  viber_df_k = list("successes"=viber_input$successes %>% filter(labels==cluster) %>% dplyr::select(-labels),
                    "trials"=viber_input$trials %>% filter(labels==cluster) %>% dplyr::select(-labels),
                    "vaf_df"=viber_input$vaf_df %>% filter(labels==cluster))
  k = viber_df_k$successes %>% nrow
  fit = ""
  try(expr = {
    fit = VIBER::variational_fit(viber_df_k$successes, viber_df_k$trials, K=k)
    fit = VIBER::choose_clusters(fit, binomial_cutoff=0, dimensions_cutoff=0, pi_cutoff=0.01)

    labels = fit$labels$cluster.Binomial
    viber_df_k$vaf_df$labels_viber = labels
    viber_df_k$vaf_df$pi_viber = fit$pi_k[labels] %>% as.vector() }, silent = T)

  try(expr = {if (fit == "") { viber_df_k$vaf_df$labels_viber = ""; viber_df_k$vaf_df$pi_viber = 0} },
      silent = T)

  return(list("df"=viber_df_k$vaf_df, "fit"=fit))
}


