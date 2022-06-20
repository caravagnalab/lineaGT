# TODO add reference and export function
# `label` is to give a label to the run
# `lineages` is to select on which lineages to run the inference,
# otherwise it will be performed on all of them
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

    print(cluster)
    # print(x.muts.k)

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

    if (x.muts.k$K * 0.01 < 1) pi_cutoff = 0.005 else pi_cutoff = 0.01
    x.muts.k = VIBER::choose_clusters(x.muts.k,
                                    binomial_cutoff=0,
                                    dimensions_cutoff=0,
                                    pi_cutoff=pi_cutoff,
                                    re_assign=T)

    labels = x.muts.k$labels$cluster.Binomial

    input.k$vaf.df$labels_viber = labels
    input.k$vaf.df$pi_viber = x.muts.k$pi_k[labels] %>% as.vector()

    if (infer_phylo)
      tree = fit_trees(x.muts.k)

  }, silent = T)

  if (purrr::is_empty(x.muts.k)) {
    input.k$vaf.df$labels_viber = "C1"
    input.k$vaf.df$pi_viber = NA
  }

  input.k$vaf.df = input.k$vaf.df %>%
    wide_to_long_muts() %>%
    add_theta_to_vaf(x.muts.k=x.muts.k, cluster=cluster)

  return(list("df"=input.k$vaf.df, "fit"=x.muts.k, "tree"=tree))
}


