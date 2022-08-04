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

fit_phylogenies = function(x, vaf.df=NULL, min_frac=0, highlight=list(), fit_muts=FALSE) {

  if (!have_muts_fit(x) && !have_vaf_df(x) && is.null(vaf.df))
    return(cli::cli_alert_warning("An input VAF dataframe is required"))

  if (!is.null(vaf.df) && have_vaf_df(x))
    cli::cli_alert_warning("Using the input mutation data but a VAF dataframe is already present in the object.")
  else if (is.null(vaf.df) && have_vaf_df(x))
    vaf.df = x %>% get_vaf_dataframe()

  clusters_joined = get_highlight(x, min_frac, highlight)
  trees = list()
  if (!have_muts_fit(x) || fit_muts)
    return(
      x %>% fit_mutations(vaf.df=vaf.df,
                          highlight=clusters_joined,
                          lineages=lineages,
                          infer_phylo=TRUE)
    )

  x.muts.all = x %>% get_muts_fit()

  if (purrr::is_empty(x.muts.all))
    return( cli::cli_alert_warning("No mutations clustering has been performed yet.
                                    Run again the function with {.field {'fit_muts'} = TRUE}") )

  trees = list()
  for (cluster in clusters_joined) {
    x.muts.k = x.muts.all[[cluster]]
    if (!is.null(x.muts.k) && x.muts.k$K > 1) {
      tt = fit_trees(x.muts.k, cluster)
      trees[[cluster]] = tt
    }
  }

  x = add_phylo(x, trees, label=label)

  return(x)
}


# to infer the tree on a single cluster
fit_trees = function(fit_viber, cluster) {
  tree = list()
  if (fit_viber$K > 1)
    tryCatch(
      expr = { tree = suppressWarnings(withTimeout(run_ctree(fit_viber, cluster), timeout=300, onTimeout="error")) },
      error = function(e) {}
    )

  return(tree)
}


run_ctree = function(viber_run, clonal) {
  if (all(is.null(viber_run$data)))
    stop("Your input object should have a data field; recreate the VIBER input.")

  if (!all(c("driver", "gene") %in% colnames(viber_run$data)))
    stop("Your data should have a logical 'driver' and 'gene' column to annotate driver events, cannot build a ctree otherwise.")

  stopifnot(inherits(viber_run, "vb_bmm"))

  patientID = ifelse(is.null(viber_run$description), clonal, viber_run$description)
  patientID = gsub(pattern = " ", replacement = ".", patientID)

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
  drivers_table = dplyr::bind_cols(drivers_table, vaf_table[which(viber_run$data$driver), , drop = F])
  drivers_table$is.clonal = FALSE
  drivers_table$is.clonal[which(drivers_table$cluster == cluster_table %>%
                                  dplyr::filter(is.clonal) %>% dplyr::pull(cluster))] = TRUE
  drivers_table = drivers_table %>% dplyr::select(patientID,
                                                  variantID, is.driver, is.clonal, cluster, colnames(cx),
                                                  dplyr::everything())

  cli::cli_process_start(paste0("Starting phylogeny inference of clone ", clonal))

  tt = ctree::ctrees(CCF_clusters=cluster_table,
                  drivers=drivers_table,
                  samples=colnames(cx),
                  patient=patientID,
                  sspace.cutoff=100,
                  store.max=50,
                  n.sampling=100)

  cli::cli_process_done()

  return(tt)
}




