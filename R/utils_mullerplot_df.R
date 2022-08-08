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

  means = x %>%
    get_mean_long() %>%
    dplyr::filter(labels %in% highlight) %>%
    dplyr::rename(pop=mean_cov) %>%
    format_means_df() %>%
    dplyr::group_by(Generation, Lineage) %>%
    dplyr::mutate(Frequency=Population/sum(Population)) %>%
    dplyr::ungroup()

  edges = get_muller_edges(x, mutations=mutations, highlight=highlight)
  if (mutations) means = get_pop_muts(x, means=means, edges=edges)

  pop_df = means %>%
    add_parent(x=x) %>%  # add common parent "P" data
    add_time_0(x=x, value=value, force=add_t0) %>%  # add initial timepoint 0

    dplyr::full_join(edges, by="Identity") %>%

    dplyr::mutate(Frequency=replace(Frequency, Generation==value & !is.na(Parent), 0)) %>%
    dplyr::mutate(Frequency=replace(Frequency, Identity=="P" & Generation==value, 1)) %>%
    dplyr::mutate(Frequency=replace(Frequency, Identity=="P" & Generation!=value, 0)) %>%

    convert_tp(timepoints_to_int=timepoints_to_int) %>%  # convert timepoints to numeric values

    dplyr::select(Identity, Generation, Lineage, Population, Frequency) %>%  #, dplyr::starts_with("lm"))
    dplyr::arrange(Identity, Lineage, Generation)

  return(pop_df)
}


get_pop_muts = function(x, means, edges) {

  means.par = means %>%
    check_fracs(edges=edges, x=x) %>%  # function to check the fractions of subclones, if linear evolutions
    dplyr::rename(Parent=Identity, Pop.par=Population, Freq.par=Frequency)

  # the means dataframe must contain also the subclones
  pop = x %>%
    get_vaf_dataframe() %>%

    dplyr::select(dplyr::starts_with("theta"), labels_mut, timepoints, lineage) %>%
    unique() %>%
    dplyr::rename(Identity=labels_mut, Lineage=lineage, Generation=timepoints) %>%
    dplyr::inner_join(edges, by="Identity") %>%
    dplyr::inner_join(means.par, by=c("Parent","Generation","Lineage")) %>%
    dplyr::mutate(Population=theta_binom*Pop.par, Frequency=theta_binom*Freq.par) %>%

    dplyr::arrange(Generation, Lineage) %>%

    # check if pop of siblings are > pop of the parent
    dplyr::group_by(Lineage, Generation, Parent) %>%
    dplyr::mutate(Pop.sum=sum(Population), Freq.sum=sum(Frequency)) %>%  # sum of the subclones populations per cluster of IS, tp, lin
    dplyr::mutate(Population=replace(Population, Pop.sum > Pop.par, Population/Pop.par),
                  Frequency=replace(Frequency, Freq.sum > Freq.par, Frequency/Freq.par)) %>%
    dplyr::ungroup() %>%

    dplyr::add_row( means %>% mutate(Parent="P") ) %>%
    dplyr::select(-dplyr::ends_with(".par"), -dplyr::ends_with(".sum")) %>%

    substract_subclonal_fracs(edges=edges) %>%

    dplyr::select(-Parent)

  return(pop)
}


check_fracs = function(means, edges, x) {
  # not clonal parents
  not_clonal = edges %>% dplyr::filter(Parent!="P") %>%
    dplyr::filter(Parent %in% Identity) %>% dplyr::pull(Parent) %>%
    unique()
  if (length(not_clonal) == 0) return(means)

  frac.par = lapply(not_clonal, get_frac_parent, means=means, edges=edges, x=x) %>%
    setNames(nm=not_clonal) %>%
    tibble::as_tibble_col() %>%  # tibble of tibbles
    tidyr::unnest(cols=c(value))  # unnest the values

  return(
    means %>%
      dplyr::add_row( frac.par )
  )
}


substract_subclonal_fracs = function(pop, edges) {
  clonal = edges %>% dplyr::filter(Parent=="P") %>%
    dplyr::pull(Identity)

  subcl.fracs = lapply(clonal, get_frac_desc, pop=pop) %>%
    setNames(nm=clonal) %>%
    tibble::as_tibble_col() %>%  # tibble of tibbles
    tidyr::unnest(cols=c(value))  # unnest the values

  return(
    pop %>%
      dplyr::full_join(subcl.fracs, by=c("Generation","Lineage","Identity")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Frequency=replace(Frequency, !is.na(Freq.sum), Frequency - Freq.sum),
                    Population=replace(Population, !is.na(Pop.sum), Population - Pop.sum)) %>%
      dplyr::ungroup()
  )
}


get_frac_parent = function(node, means, edges, x) {
  n.par = get_parent(edges, node)  # node parent

  # case 1 -> node is clonal, just return the fracs from "means"
  if (n.par == "P")
    return( means %>% dplyr::filter(Identity==node) )

  # case 2 -> node is a subclone, compute the fracs as fracs_par * theta_s
  p.frac = get_frac_parent(n.par, means, edges, x) %>%   # get the frac of the parent
    dplyr::rename(Parent=Identity, Pop.par=Population, Freq.par=Frequency)

  n.frac = x %>%
    get_vaf_dataframe() %>%
    dplyr::filter(labels_mut==node) %>%
    dplyr::select(labels_mut, timepoints, lineage, theta_binom) %>%
    unique() %>%
    dplyr::rename(Identity=labels_mut, Generation=timepoints, Lineage=lineage) %>%
    dplyr::inner_join(p.frac, by=c("Generation","Lineage")) %>%
    dplyr::mutate(Frequency=Freq.par*theta_binom, Population=Pop.par*theta_binom) %>%
    dplyr::select(Identity, Generation, Lineage, Frequency, Population)

  return(n.frac)
}


get_frac_desc = function(node, pop) {

  subcl.fracs = pop %>%
    dplyr::filter(Parent==node) %>%
    dplyr::group_by(Lineage, Generation) %>%
    dplyr::summarise(Freq.sum=sum(Frequency), Pop.sum=sum(Population),
                     Identity=Parent, .groups="keep") %>% unique() %>%
    dplyr::ungroup()

  return(subcl.fracs)
}

