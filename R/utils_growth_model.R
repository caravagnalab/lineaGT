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
                             rates.log=NULL,
                             posterior_samples.exp=NULL,
                             posterior_samples.log=NULL) {

  if (!is.null(rates.exp)) params.exp = get_growth_rates_exp(rates.exp, lineages, cluster, timepoints_to_int, posterior_samples.exp)
  if (!is.null(rates.log)) params.log = get_growth_rates_log(rates.log, lineages, cluster, timepoints_to_int, posterior_samples.log)

  if (is.null(rates.exp)) return(params.log)
  if (is.null(rates.log)) return(params.exp)

  df = dplyr::inner_join(params.exp, params.log, by=c("Lineage", "Identity"))
  return( df %>% get_growth_rates_long() )
}


posterior_samples_to_df = function(numpy_array, lineages) {
  # numpy array of dimension N_iters x L
  # returns a dataframe with nrow = L and one column with a list of samples

  if (is.null(numpy_array))
    return(data.frame(Lineage=lineages) %>% dplyr::mutate(param_samples=NA))

  numpy_array = numpy_array %>% as.matrix() %>% t() %>% as.data.frame()
  rownames(numpy_array) = lineages
  result = numpy_array %>% tibble::rownames_to_column(var="Lineage") %>%
    tidyr::pivot_longer(cols=-"Lineage") %>%
    dplyr::group_by(Lineage) %>%
    dplyr::summarise(param_samples=list(value %>% setNames(name)))
}


get_growth_rates_exp = function(rates.df, lineages, cluster, timepoints_to_int,
                                posterior_samples) {

  posterior_df = data.frame()
  if (!is.null(posterior_samples)) {
    posterior_df = posterior_samples$fitness %>% posterior_samples_to_df(lineages) %>%
      dplyr::rename(post_fitness=param_samples) %>%
      dplyr::left_join(
        posterior_samples$init_time %>% posterior_samples_to_df(lineages) %>%
          dplyr::rename(post_init_time=param_samples),
        by="Lineage") %>%
      dplyr::group_by(Lineage) %>%
      tidyr::nest(posterior_samples.exp=c(post_fitness, post_init_time))
  }

  sigma = rates.df$sigma %>% as.matrix() %>% t() %>% as.data.frame()
  colnames(sigma) = c(0, timepoints_to_int)
  rownames(sigma) = lineages

  sigma = sigma %>% tibble::rownames_to_column(var="Lineage") %>%
    tidyr::pivot_longer(cols=-"Lineage") %>%
    dplyr::group_by(Lineage) %>%
    dplyr::summarise(sigma.exp=list(value %>% setNames(name)))

  pars = data.frame(
    Lineage=lineages,
    # p_rate.exp=rates.df$parent_rate,
    fitness.exp=rates.df$fitness %>% as.numeric(),
    init_t.exp=rates.df$init_time %>% as.integer()
  ) %>%
    dplyr::mutate(rate.exp=NA,
                  p_rate.exp=ifelse(is.null(rates.df$parent_rate), NA, rates.df$parent_rate)) %>%
    dplyr::inner_join(sigma, by="Lineage")

  try(expr = {
    pars = pars %>% dplyr::mutate(p_rate.exp=as.numeric(p_rate.exp))
  })

  return(
    pars %>%

      # compute the rates for subclones
      dplyr::mutate(rate.exp=replace( rate.exp, !is.na(p_rate.exp), p_rate.exp * (1+fitness.exp) ),
                    rate.exp=replace( rate.exp, is.na(p_rate.exp), fitness.exp ),
                    # sigma.exp=list( setNames(object=rates.df$sigma, nm=c(0, unlist(timepoints_to_int))) ),
                    Identity=cluster) %>%

      dplyr::left_join(posterior_df) %>%

      tibble::as_tibble()
  )
}


get_growth_rates_log = function(rates.df, lineages, cluster, timepoints_to_int,
                                posterior_samples) {

  posterior_df = data.frame()
  if (!is.null(posterior_samples)) {
    posterior_df = posterior_samples$fitness %>% posterior_samples_to_df(lineages) %>%
      dplyr::rename(post_fitness=param_samples) %>%
      dplyr::left_join(
        posterior_samples$init_time %>% posterior_samples_to_df(lineages) %>%
          dplyr::rename(post_init_time=param_samples),
        by="Lineage") %>%
      dplyr::left_join(
        posterior_samples$carr_capac %>% posterior_samples_to_df(lineages) %>%
          dplyr::rename(post_carr_capacity=param_samples),
        by="Lineage") %>%
      dplyr::group_by(Lineage) %>%
      tidyr::nest(posterior_samples.log=c(post_fitness, post_init_time, post_carr_capacity))
  }

  sigma = rates.df$sigma %>% as.matrix() %>% t() %>% as.data.frame()
  colnames(sigma) = c(0, timepoints_to_int)
  rownames(sigma) = lineages

  sigma = sigma %>% tibble::rownames_to_column(var="Lineage") %>%
    tidyr::pivot_longer(cols=-"Lineage") %>%
    dplyr::group_by(Lineage) %>%
    dplyr::summarise(sigma.log=list(value %>% setNames(name)))

  pars = data.frame(
    Lineage=lineages,
    fitness.log=rates.df$fitness %>% as.numeric(),
    K.log=rates.df$carr_capac %>% as.integer(),
    init_t.log=rates.df$init_time %>% as.integer()
  ) %>%
    dplyr::mutate(rate.log=NA,
                  p_rate.log=ifelse(is.null(rates.df$parent_rate), NA, rates.df$parent_rate)) %>%
    dplyr::inner_join(sigma, by="Lineage")

  try(expr = {
    pars = pars %>% dplyr::mutate(p_rate.log=as.numeric(p_rate.log))
  })

  return(
    pars %>%

      dplyr::mutate(rate.log=replace( rate.log, !is.na(p_rate.log), p_rate.log * (1+fitness.log) ),
                    rate.log=replace( rate.log, is.na(p_rate.log), fitness.log ),
                    Identity=cluster
      ) %>%

      dplyr::left_join(posterior_df) %>%

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
