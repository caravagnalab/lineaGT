sort_clusters_edges = function(clusters, parents) {
  # clusters is a list of clusters
  # parents is the edges tibble
  nn = c(parents$Parent, parents$Identity) %>% setdiff("P") %>% unique() %>% length()
  nodes = c(parents$Parent, parents$Identity) %>% setdiff("P") %>% unique()

  root = setdiff(parents$Parent, parents$Identity)
  node.name = parents %>% dplyr::filter(Parent==root) %>% dplyr::pull(Identity)
  cls.sort = c(root)

  for (i in 1:nn) {
    if (length(node.name)>1) { cls.sort = c(cls.sort, node.name); break }
    if (node.name %in% clusters) cls.sort = c(cls.sort, node.name)
    node.name = parents %>% dplyr::filter(Parent==node.name) %>% dplyr::pull(Identity)
  }

  return(
    intersect(cls.sort, clusters) %>% setdiff("P") %>% unique()
  )
}


get_parent_rate = function(parents, rates.df, cluster) {

  if (!cluster %in% parents$Identity) return(list("exp"=NULL, "log"=NULL))
  if (dplyr::filter(parents, Identity==cluster)$Parent=="P") return(list("exp"=NULL, "log"=NULL))

  par = parents %>% filter(Identity==cluster) %>% dplyr::pull(Parent)  ## parent of subcl
  par.rates = rates.df %>%
    filter(Identity==par) %>%
    dplyr::select(dplyr::contains("rate.exp"), dplyr::contains("rate.log"), -dplyr::contains("p_rate")) %>%
    as.list()

  if (NA %in% unlist(par.rates)) return(list("exp"=NULL, "log"=NULL))

  return(par.rates)
}


get_growth_params = function(timepoints_to_int,
                             lineages,
                             cluster,
                             rates.exp=NULL,
                             rates.log=NULL) {

  if (!is.null(rates.exp)) params.exp = get_growth_rates_exp(rates.exp, lineages, cluster, timepoints_to_int)
  if (!is.null(rates.log)) params.log = get_growth_rates_log(rates.log, lineages, cluster, timepoints_to_int)

  if (is.null(rates.exp)) return(params.log)
  if (is.null(rates.log)) return(params.exp)
  return( dplyr::inner_join(params.exp, params.log, by=c("Lineage", "Identity")) )
}


get_growth_rates_exp = function(rates.df, lineages, cluster, timepoints_to_int) {
  return(
    data.frame() %>%

      # initialize the columns
      tibble::add_column("Lineage"=as.character(NA),
                         "fitness.exp"=NA,
                         "init_t.exp"=NA,
                         "p_rate.exp"=NA,
                         "sigma.exp"=NA,
                         "rate.exp"=as.numeric(NA)) %>%

      # add the values
      tibble::add_row(Lineage=lineages,
                      p_rate.exp=rates.df$parent_rate,
                      fitness.exp=rates.df$fitness,
                      init_t.exp=rates.df$init_time) %>%

      # compute the rates for subclones
      dplyr::mutate(rate.exp=replace( rate.exp, !is.na(p_rate.exp), p_rate.exp * (1+fitness.exp) ),
                    rate.exp=replace( rate.exp, is.na(p_rate.exp), fitness.exp ),
                    sigma.exp=list( setNames(object=rates.df$sigma, nm=c(0, unlist(timepoints_to_int))) ),
                    Identity=cluster) %>%

      tibble::as_tibble()
  )
}


get_growth_rates_log = function(rates.df, lineages, cluster, timepoints_to_int) {
  return(
    data.frame() %>%
      tibble::add_column("Lineage"=as.character(NA),
                         "fitness.log"=NA,
                         "K.log"=NA,
                         "init_t.log"=NA,
                         "p_rate.log"=NA,
                         "sigma.log"=NA,
                         "rate.log"=as.numeric(NA)) %>%

      tibble::add_row(fitness.log=rates.df$fitness,
                      K.log=rates.df$carr_capac,
                      init_t.log=rates.df$init_time,
                      Lineage=lineages,
                      p_rate.log=rates.df$parent_rate) %>%

      dplyr::mutate(rate.log=replace( rate.log, !is.na(p_rate.log), p_rate.log * (1+fitness.log) ),
                    rate.log=replace( rate.log, is.na(p_rate.log), fitness.log ),
                    sigma.log=list( setNames(object=rates.df$sigma, nm=c(0, unlist(timepoints_to_int))) ),
                    Identity=cluster
      ) %>%

      tibble::as_tibble()
  )
}


