#' Functions to run VIBER for mutations clustering
#'
#' @description add
#'
#' @param obj mvnmm object.
#' @param vaf_df a VAF dataframe, with required columns mutation,dp,ref,alt of depth, reference and alternative reads
#' @return mvnmm object with added the VAF dataframe annotated with the found clusters and the VIBER fit for each cluster
#'
#' @import VIBER
#'
#' @export run_viber
#'
#' @examples
#' obj = run_viber(obj, lineaGT::homo_ltr.vaf, min_ccf=0.07)

run_viber = function(obj, vaf_df, min_ccf=0, highlight=list()) {
  viber_df = get_input_viber(vaf_df)

  if (purrr::is_empty(highlight)) highlight = viber_df$joined$labels %>% unique() %>% droplevels()
  clusters_joined = intersect(select_relevant_clusters(obj, min_ccf), highlight)

  joined = data.frame()
  fit_all = list()
  for (cluster in clusters_joined) {
    fit_k = fit_cluster(viber_df, cluster=cluster)
    joined = rbind(joined, fit_k$df)
    fit_all[[cluster]] = fit_k$fit
  }
  obj$dataframe_vaf = joined %>% as_tibble()
  obj$viber_run = fit_all
  return(obj)
}


fit_cluster = function(viber_df, cluster, by_lineage=F) {
  viber_df_k = list("successes"=viber_df$successes %>% filter(labels==cluster) %>%
                      dplyr::select(starts_with(c("early", "mid", "late"))) %>% as_tibble(),
                    "trials"=viber_df$trials %>% filter(labels==cluster) %>%
                      dplyr::select(starts_with(c("early", "mid", "late"))) %>% as_tibble(),
                    "joined"=viber_df$joined %>% filter(labels==cluster))
  k = viber_df_k$successes %>% nrow

  try(expr = {
    fit = VIBER::variational_fit(viber_df_k$successes, viber_df_k$trials, K=k)
    fit = VIBER::choose_clusters(fit, binomial_cutoff=0, dimensions_cutoff=0, pi_cutoff=0.01)

    labels = fit$labels$cluster.Binomial
    viber_df_k$joined$labels_viber = labels
    viber_df_k$joined$pi_viber = fit$pi_k[labels] %>% as.vector() }, silent = T)

  if (!"pi_viber" %in% (viber_df_k$joined %>% colnames)) {
    fit = ""
    viber_df_k$joined$labels_viber = ""
    viber_df_k$joined$pi_viber = 0
  }

  return(list("df"=viber_df_k$joined, "fit"=fit))
}


