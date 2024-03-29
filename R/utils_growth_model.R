sort_clusters_edges = function(clusters, edges, clonal.cl=NULL) {
  if (!is.null(clonal.cl)) edges = edges %>% filter_edges(cluster=clonal.cl)

  # clusters is a list of clusters
  # edges is the edges tibble
  nn = c(edges$Parent, edges$Identity) %>% setdiff("P") %>% unique() %>% length() # the n of nodes
  nodes = c(edges$Parent, edges$Identity) %>% setdiff("P") %>% unique() # all the nodes in edges

  root = setdiff(edges$Parent, edges$Identity) # the root
  node.name = edges %>% dplyr::filter(Parent==root) %>% dplyr::pull(Identity) # the root children
  cls.sort = c(root)

  for (i in 1:nn) {

    for (ii in node.name) if (ii %in% clusters) cls.sort = c(cls.sort, ii)
    # take the lower level of the tree
    node.name = edges %>% dplyr::filter(Parent %in% node.name) %>% dplyr::pull(Identity)
  }

  return(
    intersect(cls.sort, clusters) %>% setdiff("P") %>% unique()
  )
}


filter_edges = function(edges, cluster) {
  return(
    edges %>%
      dplyr::filter(Identity==cluster |
                      Parent==cluster |
                      grepl(paste0(cluster, "."), Identity, fixed=T) |
                      grepl(paste0(cluster, "."), Parent, fixed=T))
  )
}


get_parent_rate = function(par, rates.df, cluster) {

  if (is.na(par) || par == "P") return(list("exp"=NULL, "log"=NULL))

  # if (!cluster %in% parents$Identity) return(list("exp"=NULL, "log"=NULL))
  # if (dplyr::filter(parents, Identity==cluster)$Parent=="P") return(list("exp"=NULL, "log"=NULL))

  # par = get_parent(parents, cluster)
  par.rates = rates.df %>%
    filter(Identity==par) %>%
    dplyr::select("rate", "type")

  p_rates = list("exp"=dplyr::filter(par.rates, type=="exp")$rate,
                 "log"=dplyr::filter(par.rates, type=="log")$rate)

  if (NA %in% unlist(par.rates)) return(list("exp"=NULL, "log"=NULL))

  return(p_rates)
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

  df = dplyr::inner_join(params.exp, params.log, by=c("Lineage", "Identity"))
  return( df %>% get_growth_rates_long() )
}


get_growth_rates_exp = function(rates.df, lineages, cluster, timepoints_to_int) {
  pars = data.frame() %>%

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
                    fitness.exp=rates.df$fitness %>% as.numeric(),
                    init_t.exp=rates.df$init_time %>% as.integer())

  try(expr = {
    pars = pars %>% dplyr::mutate(p_rate.exp=as.numeric(p_rate.exp))
  })

  return(
    pars %>%

      # compute the rates for subclones
      dplyr::mutate(rate.exp=replace( rate.exp, !is.na(p_rate.exp), p_rate.exp * (1+fitness.exp) ),
                    rate.exp=replace( rate.exp, is.na(p_rate.exp), fitness.exp ),
                    sigma.exp=list( setNames(object=rates.df$sigma, nm=c(0, unlist(timepoints_to_int))) ),
                    Identity=cluster) %>%

      tibble::as_tibble()
  )
}


get_growth_rates_log = function(rates.df, lineages, cluster, timepoints_to_int) {
  pars = data.frame() %>%
    tibble::add_column("Lineage"=as.character(NA),
                       "fitness.log"=NA,
                       "K.log"=NA,
                       "init_t.log"=NA,
                       "p_rate.log"=NA,
                       "sigma.log"=NA,
                       "rate.log"=as.numeric(NA)) %>%

    tibble::add_row(fitness.log=rates.df$fitness %>% as.numeric(),
                    K.log=rates.df$carr_capac %>% as.integer(),
                    init_t.log=rates.df$init_time %>% as.integer(),
                    Lineage=lineages,
                    p_rate.log=rates.df$parent_rate)

  try(expr = {
    pars = pars %>% dplyr::mutate(p_rate.log=as.numeric(p_rate.log))
  })

  return(
    pars %>%

      dplyr::mutate(rate.log=replace( rate.log, !is.na(p_rate.log), p_rate.log * (1+fitness.log) ),
                    rate.log=replace( rate.log, is.na(p_rate.log), fitness.log ),
                    sigma.log=list( setNames(object=rates.df$sigma, nm=c(0, unlist(timepoints_to_int))) ),
                    Identity=cluster
      ) %>%

      tibble::as_tibble()
  )
}


get_growth_rates_long = function(rates) {
  # rates = x %>% get_growth_rates()
  rates.exp = rates %>% dplyr::select(Lineage, Identity, dplyr::ends_with(".exp")) %>%
    dplyr::rename_with(stringr::str_replace_all, pattern=".exp", replacement="") %>%
    dplyr::mutate(K=NA, type="exp")

  rates.log = rates %>% dplyr::select(Lineage, Identity, dplyr::ends_with(".log")) %>%
    dplyr::rename_with(stringr::str_replace_all, pattern=".log", replacement="") %>%
    dplyr::mutate(type="log")

  return(
    rates.exp %>%
      dplyr::add_row(rates.log)
  )
}
