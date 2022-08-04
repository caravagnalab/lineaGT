get_muller_edges = function(x,
                            mutations=FALSE,
                            tree_score=1,
                            highlight=c()) {

  edges = data.frame("Parent"="P", "Identity"=get_unique_labels(x))
  if (!mutations) return(edges)

  edges_mut = x %>% get_parents(tree_score=tree_score)
  missed = setdiff(get_unique_muts_labels(x), edges_mut %>% dplyr::pull(Identity))

  return(edges_mut %>% dplyr::add_row(edges))
}


get_parents = function(x,
                       highlight=c(),
                       tree_score=1) {

  highlight = get_highlight(x, highlight, min_frac=0, mutations=F)

  # create an empty dataset with colnames
  edges = setNames(data.frame(matrix(ncol=2, nrow=0)), c("Parent", "Identity")) %>%
    tibble::as_tibble() %>%
    mutate(Parent=as.character(Parent), Identity=as.character(Identity))

  trees = get_trees(x)
  for (cluster in highlight) {
    tree = trees[[cluster]]
    muts = x %>% get_unique_muts_labels(cluster=cluster)

    if (!purrr::is_empty(tree))
      edges = edges %>%
        dplyr::add_row(
          tree[[tree_score]] %>%
            get_adj() %>%
            as.data.frame() %>%
            rownames_to_column(var="Label") %>%
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
                          force=T,
                          tree_score=1) {

  timepoints_to_int = map_timepoints_int(x, timepoints_to_int=timepoints_to_int)

  means = x %>%
    get_mean_long() %>%
    dplyr::rename(pop=mean_cov)

  if (mutations) {
    edges = get_muller_edges(x, highlight=highlight, tree_score=tree_score, mutations=T) %>%
      dplyr::rename(labels_mut=Identity, labels=Parent)
    # the means dataframe must contain also the subclones
    means = x %>%
      get_vaf_dataframe() %>%

      dplyr::select(dplyr::starts_with("theta"), dplyr::starts_with("labels"), timepoints, lineage) %>%
      unique() %>%
      dplyr::inner_join(edges, by=c("labels_mut","labels")) %>%
      dplyr::inner_join(means, by=c("labels", "timepoints", "lineage")) %>%
      dplyr::mutate(parent_pop=pop) %>%
      dplyr::mutate(pop=theta_binom*pop) %>%

      dplyr::rename(parent=labels, labels=labels_mut) %>%

      # check if pop of siblings are > pop of the parent
      dplyr::group_by(lineage, timepoints, parent) %>%
      dplyr::mutate(pop_sum=sum(pop)) %>%
      dplyr::mutate(pop=replace(pop, pop_sum > parent_pop, pop/pop_sum)) %>%
      dplyr::ungroup() %>%

      dplyr::add_row( means %>% mutate(parent="P") )

  }

  value = "init"
  if (!0 %in% timepoints_to_int && force)
    timepoints_to_int = c(0, timepoints_to_int) %>% setNames(nm=c(value, names(timepoints_to_int)))
  else
      value = which(timepoints_to_int==0) %>% names()

  pop_df = means %>%
    format_means_df() %>%  # to format the dataframe with correct colnames ecc
    add_parent(x=x) %>%  # add common parent "P" data
    add_time_0(x=x, value=value, force=force) %>%  # add initial timepoint 0
    dplyr::group_by(Generation, Lineage) %>%
    dplyr::mutate(Frequency=Population/sum(Population)) %>%
    dplyr::ungroup() %>%
    convert_tp(timepoints_to_int=timepoints_to_int) %>%  # convert timepoints to numeric values
    # add_exp_fit_coeff(x=x, add_exp_coef=exp_coef) %>%
    dplyr::select(Identity, Generation, Lineage, Population, Frequency) %>%  #, dplyr::starts_with("lm"))
    dplyr::arrange(Identity, Lineage, Generation)

  return(pop_df)
}


compute_mean_small_populations = function(mean, x) {

}


format_means_df = function(mean_df) {
  return(
    mean_df %>%
      dplyr::rename(Identity=labels, Generation=timepoints, Lineage=lineage, Population=pop) %>%
      dplyr::mutate(Identity=as.character(Identity)) %>%
      dplyr::mutate(Population=ifelse(Population==0, 0.001, Population))
  )
}


add_parent = function(pop_df, x=x) {
  return(
    pop_df %>%
      dplyr::add_row(
        Identity=rep( "P", times = x$data.shape[2] ),
        Population=rep( 0, times = x$data.shape[2] ),
        Generation=rep( get_timepoints(x), times = get_lineages(x) %>% length() ),
        Lineage=rep( x$lineages, each = get_timepoints(x) %>% length() )
      )
  )
}


add_time_0 = function(pop_df,
                      x=x,
                      force=T,
                      value="init") {

  if (!force) return(pop_df)  # when we do not add time 0 -> like when used to compute the clusters abundance

  n_tp = x %>% get_timepoints() %>% length()

  if (value %in% (pop_df$Generation %>% unique())) return(pop_df)

  ids = pop_df %>% filter(Identity!="P") %>% dplyr::pull(Identity) %>% unique()
  n_ids = length(ids)
  n_lins = x %>% get_lineages() %>% length()
  return(
    pop_df %>%
      dplyr::add_row(
        Identity=rep( ids, times = n_lins ),
        Population=rep( 0, times = n_ids * n_lins ),
        Generation=rep( value, times = n_ids * n_lins ),
        Lineage=rep( x$lineages, each = n_ids )
      ) %>%
      dplyr::add_row(
        Identity=rep( "P", times = n_lins ),
        Population=rep( 1, times = n_lins ),
        Generation=rep( value, times = n_lins ),
        Lineage=x$lineages
      )
  )
}


convert_tp = function(pop_df,
                      timepoints_to_int) {
  # if (is.null(timepoints_to_int)) return(pop_df)
  return(
    pop_df %>%
      mutate(Generation=timepoints_to_int[Generation])# %>%
      # mutate(Generation=unlist(Generation)) %>%
      # mutate(Generation=as.numeric(Generation))
  )
}


# add_exp_fit_coeff = function(pop_df,
#                              x,
#                              add_exp_coef=TRUE) {
#
#   if (!add_exp_coef) return(pop_df)
#   if ((pop_df$Generation %>% unique() %>% length()) == 1) return(pop_df)
#
#   generation_list = pop_df$Generation %>% unique()
#
#   return(
#     pop_df %>%
#       add_time_0(x=x, force=T, value=as.numeric(0)) %>%
#       dplyr::mutate(Generation=ifelse(Generation=="init", 0, Generation)) %>%
#       dplyr::group_by(Identity, Lineage) %>%
#       dplyr::arrange(Generation, .by_group=T) %>%
#       dplyr::mutate(lm_a=coef(lm(log1p(Population)~Generation))[1],
#                     lm_r=coef(lm(log1p(Population)~Generation))[2]) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(Generation %in% generation_list)
#   )
# }


filter_muller_df = function(df,
                            highlight=highlight) {

  return(
    df %>%
      dplyr::mutate(labels=Identity) %>%
      tidyr::separate(labels, into=c("labels", "labels_mut"), sep="[.]", fill="right") %>%
      dplyr::filter(labels %in% c(highlight, "P")) %>%
      dplyr::select(-"labels", -"labels_mut")
  )
}


pop_df_add_empty = function(mullerdf) {

  totals = mullerdf %>%
    dplyr::group_by(Generation, Lineage) %>%
    dplyr::summarise(tot=sum(Population)) %>%
    ungroup()
  max_tot = max(totals$tot)
  mullerdf$Group_id = as.character(mullerdf$Group_id)

  new_rows1 = mullerdf %>%
    dplyr::group_by(Generation, Lineage) %>%
    dplyr::summarise(Identity=NA,
                     Population=-sum(Population)/2 + 1.1 * max_tot/2) %>%
    unique() %>%
    mutate(Frequency=NA,
           Group_id="___special_empty",
           Unique_id=paste0("___special_empty_", Generation))

  new_rows2 = mullerdf %>%
    dplyr::group_by(Generation, Lineage) %>%
    dplyr::summarise(Identity=NA,
                     Population=dplyr::first(Population)) %>%
    mutate(Frequency=NA,
           Group_id="___special_emptya",
           Unique_id=paste0("___special_emptya_",Generation))

  mullerdf = new_rows1 %>% dplyr::bind_rows(mullerdf) %>%
    dplyr::arrange(Generation) %>%
    ungroup() %>%
    dplyr::bind_rows(new_rows2) %>%
    dplyr::arrange(Generation) %>%
    ungroup() %>%
    dplyr::group_by(Generation, Lineage) %>%
    dplyr::mutate(Frequency=Population/sum(Population)) %>%
    ungroup()

  mullerdf$Group_id = factor(mullerdf$Group_id,
                             levels=rev(unlist(
                               as.data.frame(mullerdf %>%
                                               dplyr::filter_(~Generation == max(Generation)) %>%
                                               dplyr::select_(~Group_id)),
                               use.names=FALSE) %>% unique()))
  return(mullerdf)
}


