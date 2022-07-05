#' Fit the phylogenetic trees
#'
#' @param x add
#' @param vaf.df add
#' @param min_frac add
#' @param highlight add
#' @param do_filter add
#' @param label add
#' @param fit_viber add
#' @param lineages add
#'
#' @return
#' @export fit_phylogenies

fit_phylogenies = function(x, vaf.df=NULL, min_frac=0, highlight=list(), do_filter=FALSE,
                           label="", fit_viber=FALSE, lineages=c()) {

  if (is.null(vaf.df) && !"vaf.dataframe" %in% names(x))
    message("A dataframe with the mutations is required!")
  else if (is.null(vaf.df))
    vaf.df = x %>% get_vaf_dataframe(label=label)

  clusters_joined = get_highlight(x, min_frac, highlight)
  trees = list()
  if (!"viber_run" %in% names(x) || fit_viber)
    return(
      x %>% fit_mutations(vaf.df=vaf.df,
                      highlight=clusters_joined,
                      lineages=lineages,
                      label=label,
                      infer_phylo=TRUE)
    )

  viber_run_all = x %>% get_muts_fit(label=label)

  if (is.null(viber_run_all))
    return(message("No mutations clustering has been performed with the input label!"))

  for (cluster in clusters_joined) {
    viber_run = viber_run_all[[cluster]]
    tt = fit_trees(viber_run, cluster)
    trees[[cluster]] = tt
  }

  x = add_phylo(x, trees, label=label)

  return(x)
}


# to infer the tree on a single cluster
fit_trees = function(fit_viber, clonal) {
  if (length(fit_viber$labels$cluster.Binomial %>% unique) > 1)
    tree = run_ctree(fit_viber, clonal)
  else
    tree = list()

  return(tree)
}


run_ctree = function(viber_run, clonal) {
  # viber_run here is a viber fit
  if (all(is.null(viber_run$data)))
    stop("Your input object should have a data field; recreate the VIBER input.")
  if (!all(c("driver", "gene") %in% colnames(viber_run$data)))
    stop("Your data should have a logical 'driver' and 'gene' column to annotate driver events, cannot build a ctree otherwise.")
  stopifnot(inherits(viber_run, "vb_bmm"))
  patientID = ifelse(is.null(viber_run$description), "VIBER dataset",
                     viber_run$description)
  patientID = gsub(pattern = " ", replacement = "_", patientID)
  pi = viber_run$pi_k[((viber_run$N * viber_run$pi_k) %>% round) > 0]
  theta = viber_run$theta_k[, names(pi), drop = T]
  cluster_table = data.frame(cluster = colnames(theta), stringsAsFactors = FALSE) %>%
    tibble::as_tibble()
  cluster_table = dplyr::bind_cols(cluster_table, t(theta) %>% tibble::as_tibble())
  cluster_table$nMuts = table(viber_run$labels)[cluster_table$cluster] %>%
    as.vector()
  clonal_cluster = apply(theta, 1, which.max)
  clonal_cluster = colnames(theta)[clonal_cluster]

  cx = viber_run$x %>% dplyr::select(-starts_with("cluster"))
  cy = viber_run$y %>% dplyr::select(-starts_with("cluster"))

  clonal_cluster = clonal
  theta = theta %>% as.data.frame()
  theta$P = 1.0
  cluster_table = cluster_table %>%
    dplyr::add_row(cluster=clonal, nMuts=1) %>%
    replace(is.na(.), 1)
  cx = cx %>% dplyr::add_row() %>% replace(is.na(.), max(cy))
  cy = cy %>% dplyr::add_row() %>% replace(is.na(.), max(cy))

  cluster_table$is.clonal = FALSE
  cluster_table$is.clonal[cluster_table$cluster %in% clonal_cluster] = TRUE
  viber_run$data$cluster = paste(unlist(viber_run$labels))
  drivers_collapse = viber_run$data %>% dplyr::filter(driver) %>% pull(cluster) %>%
    unique
  cluster_table$is.driver = FALSE
  cluster_table$is.driver[which(cluster_table$cluster %in%
                                  drivers_collapse)] = TRUE

  vaf_table = cx/cy
  drivers_table = viber_run$data %>% tibble::as_tibble() %>% dplyr::filter(driver) %>%
    dplyr::rename(variantID = gene, is.driver = driver) %>%
    dplyr::mutate(patientID = patientID)
  drivers_table = dplyr::bind_cols(drivers_table, vaf_table[which(viber_run$data$driver),
                                                     , drop = F])
  drivers_table$is.clonal = FALSE
  drivers_table$is.clonal[which(drivers_table$cluster == cluster_table %>%
                                  dplyr::filter(is.clonal) %>% dplyr::pull(cluster))] = TRUE
  drivers_table = drivers_table %>% dplyr::select(patientID,
                                                  variantID, is.driver, is.clonal, cluster, colnames(cx),
                                                  dplyr::everything())
  tt = ctree::ctrees(CCF_clusters=cluster_table,
                     drivers=drivers_table,
                     samples=colnames(cx),
                     patient=patientID)
  return(tt)
}




