# `label` is to give a label to the run
# `lineages` is to select on which lineages to run the inference,
# otherwise it will be performed on all of them

#' Fit the mutations clustering
#'
#' @param x add
#' @param vaf.df add
#' @param infer_phylo add
#' @param min_frac add
#' @param highlight add
#' @param do_filter add
#' @param lineages add
#' @param label add
#'
#' @return
#' @export fit_mutations

fit_mutations = function(x,
                         vaf.df=NULL,
                         infer_phylo=TRUE,
                         min_frac=0,
                         highlight=list(),
                         do_filter=FALSE,
                         lineages=c(),
                         label="") {

  if (is.null(vaf.df) && !"vaf.dataframe" %in% names(x))
    message("A dataframe with the mutations is required!")
  else if (is.null(vaf.df))
    vaf.df = x %>% get_vaf_dataframe(label)

  clusters_joined = get_highlight(x, min_frac, highlight)
  x = x %>%
    annotate_vaf_df(vaf.df=vaf.df, min_frac=min_frac, label=label) %>%  # add cluster to each mutation
    check_dp(label=label)  # check for too low values

  input_viber = x %>%
    get_input_viber(lineages=lineages,
                    label=label)  # get the input to run viber

  # vaf.df.fit -> long format dataframe with already `pi`, `labels`, `theta`
  vaf.df.fit = data.frame(); x.muts = list()
  if (infer_phylo)
    x.trees = list()

  for (cluster in clusters_joined) {
    x.muts.k = fit_cluster_viber(input_viber,
                                 cluster=cluster,
                                 infer_phylo=infer_phylo)

    vaf.df.fit = rbind(vaf.df.fit, x.muts.k$df)
    x.muts[[cluster]] = x.muts.k$fit
    if (infer_phylo)
      x.trees[[cluster]] = x.muts.k$tree
  }

  x = add_muts_fit(x, x.muts=x.muts, label=label)

  x = add_vaf(x, vaf.df=vaf.df.fit, label=label)

  color_palette = update_color_palette(x, clusters_joined, label=label)
  x = add_color_palette(x, color_palette, label=label)

  x = add_phylo(x, x.trees, label=label)

  return(x)
}


fit_cluster_viber = function(input, cluster, infer_phylo=TRUE) {
  input.k = input %>% filter_viber_input(cluster=cluster)

  k = input.k$successes %>% nrow  # max n of clusters
  x.muts.k = tree = list()

  try(expr = {
    data_annotations = get_data_annotation(k)
    x.muts.k = VIBER::variational_fit(input.k$successes,
                                      input.k$trials,
                                      K=k,
                                      data=data_annotations)

    # if (x.muts.k$K * 0.01 < 1) pi_cutoff = 0.005 else pi_cutoff = 0.01
    # x.muts.k = VIBER::choose_clusters(x.muts.k,
    #                                   binomial_cutoff=0,
    #                                   dimensions_cutoff=0,
    #                                   pi_cutoff=0,
    #                                   re_assign=T)

    x.muts.k = x.muts.k %>%
      replace_labels_muts(pattern="C", replacement="S")

    labels = x.muts.k$labels$cluster.Binomial

    input.k$vaf.df$labels_viber = labels
    input.k$vaf.df$pi_viber = x.muts.k$pi_k[labels] %>% as.vector()

    if (infer_phylo)
      tree = fit_trees(x.muts.k, cluster)

  }, silent = T)

  if (purrr::is_empty(x.muts.k)) {
    input.k$vaf.df$labels_viber = "S1"
    input.k$vaf.df$pi_viber = NA
  }

  input.k$vaf.df = input.k$vaf.df %>%
    wide_to_long_muts() %>%
    add_theta_to_vaf(x.muts.k=x.muts.k, cluster=cluster)

  return(list("df"=input.k$vaf.df, "fit"=x.muts.k, "tree"=tree))
}



replace_labels_muts = function(x.muts.k, pattern, replacement) {
  x.muts.k$x$cluster.Binomial = x.muts.k$x$cluster.Binomial %>%
    str_replace_all(pattern, replacement)
  x.muts.k$y$cluster.Binomial = x.muts.k$y$cluster.Binomial %>%
    str_replace_all(pattern, replacement)

  names(x.muts.k$pi_k) = names(x.muts.k$pi_k) %>%
    str_replace_all(pattern, replacement)

  colnames(x.muts.k$theta_k) = colnames(x.muts.k$theta_k) %>%
    str_replace_all(pattern, replacement)

  names(x.muts.k$alpha) = names(x.muts.k$alpha) %>%
    str_replace_all(pattern, replacement)
  names(x.muts.k$alpha_0) = names(x.muts.k$b_0) %>%
    str_replace_all(pattern, replacement)

  colnames(x.muts.k$r_nk) = colnames(x.muts.k$r_nk) %>%
    str_replace_all(pattern, replacement)

  x.muts.k$labels$cluster.Binomial = x.muts.k$labels$cluster.Binomial %>%
    str_replace_all(pattern, replacement)

  colnames(x.muts.k$a) = colnames(x.muts.k$a) %>%
    str_replace_all(pattern, replacement)
  colnames(x.muts.k$a_0) = colnames(x.muts.k$a_0) %>%
    str_replace_all(pattern, replacement)

  colnames(x.muts.k$b) = colnames(x.muts.k$b) %>%
    str_replace_all(pattern, replacement)
  colnames(x.muts.k$b_0) = colnames(x.muts.k$b_0) %>%
    str_replace_all(pattern, replacement)

  return(x.muts.k)
}

