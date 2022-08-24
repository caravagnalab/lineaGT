## Get muller edges dataframe ####

get_muller_edges = function(x,
                            mutations=FALSE,
                            tree_score=1,
                            highlight=c()) {

  highlight = get_highlight(x, highlight=highlight, mutations=F)

  edges = data.frame("Parent"="P", "Identity"=get_unique_labels(x)) %>%
    dplyr::filter(Identity %in% highlight)
  if (!mutations) return(edges)

  edges_mut = x %>% get_edges_muts(tree_score=tree_score, highlight=highlight)
  missed = setdiff(get_unique_muts_labels(x), edges_mut %>% dplyr::pull(Identity))

  return(edges_mut %>% dplyr::add_row(edges))
}


get_edges_muts = function(x, highlight=c(), tree_score=1) {

  highlight = get_highlight(x, highlight, min_frac=0, mutations=F)

  # create an empty dataset with colnames
  edges = setNames(data.frame(matrix(ncol=2, nrow=0)), c("Parent", "Identity")) %>%
    tibble::as_tibble() %>%
    mutate(Parent=as.character(Parent), Identity=as.character(Identity))

  trees = get_trees(x)
  for (cluster in highlight) {
    tree = trees[[cluster]]
    muts = x %>% get_unique_muts_labels(cluster=cluster)

    if (!is.null(tree))
      edges = edges %>%
      dplyr::add_row(
        tree[[tree_score]] %>%
          get_adj() %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var="Label") %>%
          reshape2::melt(id="Label", variable.name="Identity") %>%
          dplyr::filter(value==1) %>%
          dplyr::filter(!Label %in% c("GL"), !Identity %in% c("GL", "P", cluster)) %>%
          dplyr::select(-value) %>%
          dplyr::rename(Parent=Label) %>%
          dplyr::mutate(Parent=ifelse(Parent==cluster, Parent, paste(cluster, Parent, sep="."))) %>%
          dplyr::mutate(Identity=paste(cluster,Identity,sep=".")) %>%
          tibble::as_tibble()
      )
    else if (length(muts) > 0)
      edges = edges %>%
      dplyr::add_row(
        data.frame(Parent=cluster, Identity=muts)
      )
  }
  return(edges)
}


## Get muller population dataframe ####

# means format must be a dataframe with columns: labels, timepoints, lineage, mean_cov
get_muller_pop = function(x,
                          mutations=F,
                          timepoints_to_int=c(),
                          highlight=c(),
                          add_t0=T,
                          tree_score=1) {

  highlight = get_highlight(x, mutations=mutations, highlight=highlight)

  # if highlight is specified -> we want to observe the frequencies of them related one to the other
  timepoints_to_int = map_timepoints_int(x, timepoints_to_int=timepoints_to_int)
  value = "init"
  if (!0 %in% timepoints_to_int && add_t0)
    timepoints_to_int = c(0, timepoints_to_int) %>% setNames(nm=c(value, names(timepoints_to_int))) else
    value = which(timepoints_to_int==0) %>% names()

  edges = get_muller_edges(x, mutations=mutations, highlight=highlight)

  means.clonal = x %>%
    get_mean_long() %>%
    dplyr::filter(labels %in% highlight) %>%
    dplyr::rename(Identity=labels, Generation=timepoints, Lineage=lineage, Population=mean_cov) %>%
    dplyr::mutate(Identity=as.character(Identity)) %>%
    dplyr::mutate(Population=ifelse(Population==0, 0.001, Population)) %>%
    dplyr::arrange(Identity, Lineage, Generation) %>%
    dplyr::mutate(theta_binom=1, theta.par=NA, Parent="P")


  pop_df = x %>%
    get_pop_muts(means=means.clonal, edges=edges, mutations=mutations) %>%
    dplyr::select(-Parent) %>%
    add_parent(x=x) %>%
    add_time_0(x=x, force=add_t0, value=value) %>%
    convert_tp(timepoints_to_int=timepoints_to_int) %>%
    dplyr::full_join(edges, by="Identity")

  return(pop_df)
}


get_pop_muts = function(x, means, edges, mutations=F) {
  if (!mutations)
    return(
      means %>%
        dplyr::select(-theta.par) %>%

        dplyr::group_by(Generation, Lineage) %>%
        dplyr::mutate(Frequency=Population/sum(Population)) %>%
        dplyr::ungroup() %>%

        dplyr::mutate(Pop.plot=Population)
      )

  # obtain the correct theta for each parent
  # columns Identity, Generation, Lineage, theta_binom, theta.par, Parent
  pop = means %>%
    correct_theta(edges=edges, x=x) %>%
    add_population(means=means, x=x, edges=edges) %>%
    dplyr::inner_join(edges, by="Identity") %>%

    group_by(Generation, Lineage) %>%
    dplyr::mutate(Frequency=ifelse(sum(Population)>0, Population / sum(Population), 0)) %>%
    dplyr::ungroup()

  pop.subcl = pop %>%
    dplyr::group_by(Parent, Generation, Lineage) %>%
    dplyr::summarise(Pop.subcl=sum(Population), .groups="keep") %>%
    dplyr::ungroup() %>%
    dplyr::filter(Parent != "P") %>%
    dplyr::rename(Identity=Parent)


  return(
    pop %>%
      dplyr::full_join(pop.subcl, by=c("Identity", "Generation", "Lineage")) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(Pop.plot=Population - Pop.subcl)
  )

  return(pop)
}


correct_theta = function(x, means, edges) {
  not.clonal = get_parents(edges, clonal=F)
  clonal = get_parents(edges, clonal=T)
  not.parents = edges$Identity %>% unique()

  # theta of all clonal clusters of IS, with theta=1 and parent="P"
  theta.df = means %>%
    dplyr::select(-Population)

  for (node in clonal) {
    theta.df = correct_theta_node(x, theta.df, edges, node)
  }

  for (node in not.clonal) {
    theta.df = correct_theta_node(x, theta.df, edges, node)
  }

  return(
    theta.df %>%
      # dplyr::rename(theta_binom=theta.par)
      dplyr::select(-theta.par)
    )
}


# given "node", it checks if its children thetas sum is lower than or equal to its own theta
correct_theta_node = function(x, theta.df, edges, node) {
  # node children
  node.c = edges %>%
    dplyr::filter(Parent==node) %>%
    dplyr::pull(Identity) %>% unique()

  if (node %in% theta.df$Identity) present = TRUE else present = FALSE

  if (all(node.c %in% theta.df$Identity))
    return(theta.df)
  if (any(node.c %in% theta.df$Identity))
    cli::cli_alert_warning("Only some descendants already in the dataframe, cluster {node}")

  theta.node = get_theta_node(x, theta.df, node, edges) %>%
    dplyr::rename(theta.par=theta_binom, Parent=Identity)

  theta.tmp = x %>%
    get_vaf_dataframe() %>%
    dplyr::rename(Identity=labels_mut, Generation=timepoints, Lineage=lineage) %>%
    dplyr::select(Identity, Generation, Lineage, theta_binom) %>%
    dplyr::inner_join(edges, by="Identity") %>%
    dplyr::filter(Parent==node) %>%
    unique() %>%

    dplyr::inner_join(theta.node, by=c("Parent","Generation","Lineage")) %>%

    dplyr::group_by(Generation, Lineage) %>%
    dplyr::mutate(theta_binom=replace(theta_binom,
                                      sum(theta_binom) > theta.par,
                                      theta_binom / sum(theta_binom) * theta.par)) %>%
    dplyr::ungroup()

  if (!present)
    theta.tmp = theta.tmp %>%
      dplyr::add_row(
        theta.node %>%
          dplyr::rename(Identity=Parent, theta_binom=theta.par)
      )

  return(
    theta.df %>%
      dplyr::add_row(
        theta.tmp
      )
  )
}


get_theta_node = function(x, theta.df, node, edges) {
  if (node %in% theta.df$Identity)
    return(
      theta.df %>%
        dplyr::filter(Identity==node) %>%
        dplyr::select(theta_binom, Identity, Generation, Lineage)
    )

  return(
    x %>%
      get_vaf_dataframe() %>%
      dplyr::filter(labels_mut==node) %>%
      dplyr::select(theta_binom, labels_mut, timepoints, lineage) %>%
      unique() %>%
      dplyr::rename(Identity=labels_mut, Generation=timepoints, Lineage=lineage)
  )
}


get_parents = function(edges, clonal) {
  if (clonal)
    return(
      edges %>%
        dplyr::filter(Parent=="P") %>%
        dplyr::pull(Identity) %>%
        unique()
    )

  return(
    edges %>%
      dplyr::filter(Parent!="P") %>%
      dplyr::filter(Parent %in% Identity) %>%
      dplyr::pull(Parent) %>%
      unique()
  )
}


add_population = function(x, means, theta.df, edges) {
  parents = edges %>% dplyr::filter(Parent!="P") %>% dplyr::pull(Parent) %>% unique()

  means = means %>%
    dplyr::select(-Parent, -theta.par)

  means.par = means %>%
    dplyr::select(-theta_binom) %>%
    dplyr::rename(Parent=Identity, Pop.par=Population)

  # first compute the population size of the parents
  pop.parents = theta.df %>%
    dplyr::filter(Identity %in% parents) %>%
    dplyr::full_join(means, by=c("Identity", "Generation", "Lineage", "theta_binom")) %>%
    dplyr::left_join(means.par, by=c("Parent", "Generation", "Lineage")) %>%

    dplyr::rowwise() %>%
    dplyr::mutate(Population=replace(Population, is.na(Population), theta_binom * Pop.par)) %>%
    dplyr::ungroup() %>%

    dplyr::select(Identity, Generation, Lineage, theta_binom, Population) %>%
    dplyr::rename(Pop.par=Population, Parent=Identity, theta.par=theta_binom)

  pop = theta.df %>%
    # add to each clonal Identity its population
    dplyr::full_join(means, by=c("Identity", "Generation", "Lineage", "theta_binom")) %>%

    # add to each subclonal parent its parent population and compute its own pop
    dplyr::left_join(pop.parents, by=c("Parent", "Generation", "Lineage")) %>%

    # compute the population size
    dplyr::rowwise() %>%
    dplyr::mutate(Population=replace(Population, is.na(Population), theta_binom * Pop.par)) %>%
    dplyr::ungroup() %>%

    dplyr::select(Identity, Generation, Lineage, Population, theta_binom)

  return( pop )
}


# Checks the fraction of non-clonal parents
# -> their pop and freqs are computed as pop_s = theta_s*pop_p (same for freq)
# -> the pop and freqs of their children are computed as theta_c*pop_s
# check_fracs = function(means, edges, x) {
#   # not clonal parents
#   not_clonal = edges %>% dplyr::filter(Parent!="P") %>%
#     dplyr::filter(Parent %in% Identity) %>% dplyr::pull(Parent) %>%
#     unique()
#   if (length(not_clonal) == 0) return(means)
#
#   frac.par = lapply(not_clonal, get_frac_parent, means=means, edges=edges, x=x) %>%
#     setNames(nm=not_clonal) %>%
#     tibble::as_tibble_col() %>%  # tibble of tibbles
#     tidyr::unnest(cols=c(value))  # unnest the values
#
#   return(
#     means %>%
#       dplyr::add_row( frac.par )
#   )
# }


# Function to retrieve the pops and fracs of the parent of "node"
# recursively retrieves/computes the pops/fracs of the parents
# get_frac_parent = function(node, means, edges, x) {
#   n.par = get_parent(edges, node)  # parent of "node"
#
#   # case 1 -> node is clonal, just return the fracs from "means"
#   if (n.par == "P")
#     return( means %>% dplyr::filter(Identity==node) )
#
#   # case 2 -> node is a subclone, compute the fracs as fracs_par * theta_s
#   # recursively retrieve the parent fractions
#   p.frac = get_frac_parent(n.par, means, edges, x) %>%   # get the frac of the parent
#     dplyr::rename(Parent=Identity, Pop.par=Population, Freq.par=Frequency, theta.par=theta_binom)
#
#   n.frac = x %>%
#     get_vaf_dataframe() %>%
#     dplyr::select(labels_mut, timepoints, lineage, theta_binom) %>%
#     dplyr::filter(labels_mut==node) %>%
#     dplyr::rename(Identity=labels_mut, Generation=timepoints, Lineage=lineage) %>%
#     unique() %>%
#
#     dplyr::inner_join(p.frac, by=c("Generation","Lineage")) %>%
#
#     # correct for the thetas > than the theta of the parent
#     dplyr::mutate(theta_binom=replace(theta_binom, theta_binom > theta.par, theta.par)) %>%
#
#     # dplyr::mutate(Frequency=Freq.par*theta_binom, Population=Pop.par*theta_binom) %>%
#     dplyr::mutate(Population=Pop.par*theta_binom) %>%
#     dplyr::select(Identity, Generation, Lineage, Frequency, Population, theta_binom)
#
#   return(n.frac)
# }


# check_thetas = function(means) {
#   tmp = means %>%
#     dplyr::mutate(cond=theta_binom>theta.par)
#
#   if (any(tmp$cond))
#     cli::cli_alert_warning("Subclone frequency > parent frequency. Resetting it to the parent frequency.")
# }


# substract_subclonal_fracs = function(pop, edges) {
#   # get the list of parents
#   parents = edges %>% dplyr::filter(Parent!="P") %>%
#     dplyr::pull(Parent) %>% unique()
#
#   # for each parent, we have two columns Pop.sum and Freq.sum
#   # -> sum of fracs and pops of all descendants
#   subcl.fracs = lapply(parents, get_frac_desc, pop=pop) %>%
#     setNames(nm=parents) %>%
#     tibble::as_tibble_col() %>%  # tibble of tibbles
#     tidyr::unnest(cols=c(value))  # unnest the values
#
#   return(
#     pop %>%
#       # it will be NA if Identity is not a parent
#       dplyr::full_join(subcl.fracs, by=c("Generation","Lineage","Identity")) %>%
#       dplyr::rowwise() %>%
#       dplyr::mutate(Freq.plot=replace(Frequency, !is.na(Freq.sum) & Frequency>=Freq.sum, Freq.plot - Freq.sum),
#                     Pop.plot=replace(Population, !is.na(Pop.sum) & Population>=Pop.sum, Pop.plot - Pop.sum)) %>%
#
#       dplyr::ungroup()
#   )
# }


# get_frac_desc = function(node, pop) {
#
#   subcl.fracs = pop %>%
#     dplyr::filter(Parent==node) %>%
#     dplyr::group_by(Lineage, Generation) %>%
#     dplyr::summarise(Freq.sum=sum(Frequency), Pop.sum=sum(Population),
#                      Identity=Parent, .groups="keep") %>% unique() %>%
#     dplyr::ungroup()
#
#   return(subcl.fracs)
# }

