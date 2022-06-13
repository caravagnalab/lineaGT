get_muller_edges = function(x, mutations=FALSE) {
  edges = data.frame("Parent"="P", "Identity"=get_unique_labels(x))

  if (!mutations) return(edges)

  edges = x %>% get_binomial_theta() %>%
    inner_join(x %>% get_mean_long(), by=c("labels", "timepoints", "lineage")) %>%
    dplyr::rename(Parent=labels, Identity=labels_mut) %>%
    dplyr::select(Parent, Identity) %>% unique() %>%
    dplyr::add_row(edges)

  return(
    x %>%
      get_parents() %>%
      dplyr::full_join(edges, by=c("Parent", "Identity")) %>%
      mutate(Parent=ifelse(is.na(Label), Parent, Label)) %>%
      dplyr::select(-Label)
  )
}


get_parents = function(x, highlight=c()) {
  if (purrr::is_empty(highlight))
    highlight = get_unique_labels(x)

  edges = setNames(data.frame(matrix(ncol=3, nrow=0)),
                   c("Parent", "Identity", "Label")) %>%
    tibble::as_tibble() %>%
    mutate(Parent=as.character(Parent), Identity=as.character(Identity), Label=as.character(Label))

  for (cluster in highlight) {
    tree = x$trees[[cluster]]

    if (!purrr::is_empty(tree)) {
      edges = edges %>%
        dplyr::add_row(
          tree[[1]] %>%
            get_adj() %>%
            as.data.frame() %>%
            rownames_to_column(var="Label") %>%
            reshape2::melt(id="Label", variable.name="Identity") %>%
            filter(value==1) %>%
            filter(Label != "GL", Identity != "GL") %>%
            dplyr::select(-value) %>%
            dplyr::mutate(Label=paste(cluster, Label, sep="."),
                          Identity=paste(cluster, Identity, sep="."),
                          Parent=cluster) %>%
            tibble::as_tibble()
        )

    }
  }
  return(edges)
}


get_adj = function(tree) {
  return(
    tree$adj_mat
  )
}


# means format must be a dataframe with columns: labels, timepoints, lineage, mean_cov
get_muller_pop = function(x, map_tp_time=list("init"=0,"early"=60,"mid"=140,"late"=280),
                          mutations=FALSE) {
  means = x %>% get_mean_long()

  if (mutations)
    # the means dataframe must contain also the subclones
    means = x %>% get_binomial_theta() %>%
      inner_join(means, by=c("labels", "timepoints", "lineage")) %>%
      mutate(mean_cov=theta/100*mean_cov) %>%
      dplyr::rename(parent=labels, labels=labels_mut) %>%
      dplyr::add_row( means %>% mutate(parent="P") )

  pop_df = means %>%
    format_means_df() %>%  # to format the dataframe with correct colnames ecc
    add_parent(x=x) %>%  # add common parent "P" data
    # add_time_0() %>%
    convert_tp(mapping=map_tp_time) %>%  # convert timepoints to numeric values
    add_exp_fit_coeff() %>%
    dplyr::select(Identity, Generation, Lineage, Population, Frequency, lm_a, lm_r)

  return(pop_df)
}


format_means_df = function(mean_df) {
  return(
    mean_df %>%
      dplyr::rename(Identity=labels, Generation=timepoints, Lineage=lineage, Population=mean_cov) %>%
      dplyr::mutate(Identity=as.character(Identity)) %>%
      dplyr::mutate(Population=ifelse(Population==0, 0.001, Population)) %>%
      dplyr::group_by(Generation, Lineage) %>%
      dplyr::mutate(Frequency=Population/sum(Population)) %>%
      dplyr::ungroup()
  )
}

add_parent = function(pop_df, x=x) {
  return(
    pop_df %>%
      dplyr::add_row(
        Identity=rep( "P", times = x$`T` ),
        Population=rep( 1, times = x$`T` ),
        Frequency=rep( 1, times = x$`T` ),
        Generation=rep( get_timepoints(x), times = get_lineages(x) %>% length() ),
        Lineage=rep( x$lineages, each = get_timepoints(x) %>% length() )
      )
  )
}

add_time_0 = function(pop_df, x=x) {
  ids = pop_df$Identity %>% unique()
  n_ids = length(ids)
  return(
    pop_df %>%
      dplyr::add_row(
        Identity=rep( ids, times = get_lineages(x) %>% length() ),
        Population=rep( 0, times = n_ids * get_lineages(x) %>% length() ),
        Frequency=rep( 0, times = n_ids * get_lineages(x) %>% length() ),
        Generation=rep( "init", times = n_ids * get_lineages(x) %>% length() ),
        Lineage=rep( x$lineages, each = n_ids )
      )
  )
}

convert_tp = function(pop_df, mapping=list("init"="0","early"="60","mid"="140","late"="280")) {
  return(
    pop_df %>%
      mutate(Generation=mapping[Generation]) %>%
      mutate(Generation=unlist(Generation)) %>%
      mutate(Generation=as.numeric(Generation))
  )
}

add_exp_fit_coeff = function(pop_df) {
  return(
    pop_df %>%
      dplyr::group_by(Identity, Lineage) %>%
      dplyr::mutate(lm_a=coef(lm(log1p(Population)~Generation))[1],
                    lm_r=coef(lm(log1p(Population)~Generation))[2]) %>%
      dplyr::ungroup()
  )
}



