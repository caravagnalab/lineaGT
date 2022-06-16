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




run_ctree = function(x) {
  # x here is a viber fit
  if (all(is.null(x$data)))
    stop("Your input object should have a data field; recreate the VIBER input.")
  if (!all(c("driver", "gene") %in% colnames(x$data)))
    stop("Your data should have a logical 'driver' and 'gene' column to annotate driver events, cannot build a ctree otherwise.")
  stopifnot(inherits(x, "vb_bmm"))
  patientID = ifelse(is.null(x$description), "VIBER dataset",
                     x$description)
  patientID = gsub(pattern = " ", replacement = "_", patientID)
  pi = x$pi_k[((x$N * x$pi_k) %>% round) > 0]
  theta = x$theta_k[, names(pi), drop = T]
  cluster_table = data.frame(cluster = colnames(theta), stringsAsFactors = FALSE) %>%
    as_tibble()
  cluster_table = bind_cols(cluster_table, t(theta) %>% as_tibble)
  cluster_table$nMuts = table(x$labels)[cluster_table$cluster] %>%
    as.vector()
  clonal_cluster = apply(theta, 1, which.max)
  clonal_cluster = colnames(theta)[clonal_cluster]
  if (clonal_cluster %>% unique() %>% length() == 1)
    clonal_cluster = which.max(table(clonal_cluster)) %>% names
  else{
    clonal_cluster = "P"
    theta = theta %>% as.data.frame()
    theta$P = 1.0
  }
  cli::cli_text("Estimated clonal cluster {.value {clonal_cluster}} from VIBER fit.")
  cluster_table$is.clonal = FALSE
  cluster_table$is.clonal[cluster_table$cluster %in% clonal_cluster] = TRUE
  x$data$cluster = paste(unlist(x$labels))
  drivers_collapse = x$data %>% dplyr::filter(driver) %>% pull(cluster) %>%
    unique
  cluster_table$is.driver = FALSE
  cluster_table$is.driver[which(cluster_table$cluster %in%
                                  drivers_collapse)] = TRUE
  cli::cli_text("Found {.value {sum(cluster_table$is.driver)}} driver event(s) in VIBER fits.")
  cx = x$x %>% dplyr::select(-starts_with("cluster"))
  cy = x$y %>% dplyr::select(-starts_with("cluster"))
  vaf_table = cx/cy
  drivers_table = x$data %>% as_tibble() %>% dplyr::filter(driver) %>%
    dplyr::rename(variantID = gene, is.driver = driver) %>%
    dplyr::mutate(patientID = patientID)
  drivers_table = bind_cols(drivers_table, vaf_table[which(x$data$driver),
                                                     , drop = F])
  drivers_table$is.clonal = FALSE
  drivers_table$is.clonal[which(drivers_table$cluster == cluster_table %>%
                                  dplyr::filter(is.clonal) %>% dplyr::pull(cluster))] = TRUE
  drivers_table = drivers_table %>% dplyr::select(patientID,
                                                  variantID, is.driver, is.clonal, cluster, colnames(cx),
                                                  dplyr::everything())
  trees = ctree::ctrees(CCF_clusters = cluster_table, drivers = drivers_table,
                        samples = colnames(cx), patient = patientID)
  return(trees)
}

