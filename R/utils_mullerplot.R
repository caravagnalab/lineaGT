get_muller_edges = function(x, mutations=FALSE, label="", tree_score=1) {
  edges = data.frame("Parent"="P", "Identity"=get_unique_labels(x))

  if (!mutations) return(edges)

  edges = x %>%
    get_vaf_dataframe(label=label) %>%
    inner_join(x %>% get_mean_long(), by=c("labels", "timepoints", "lineage")) %>%
    dplyr::rename(Parent=labels, Identity=labels_mut) %>%
    dplyr::select(Parent, Identity) %>% unique() %>%
    dplyr::add_row(edges)

  return(
    x %>%
      get_parents(label=label, tree_score=tree_score) %>%
      dplyr::full_join(edges, by=c("Parent", "Identity")) %>%
      mutate(Parent=ifelse(is.na(Label), Parent, Label)) %>%
      dplyr::select(-Label)
  )
}


get_parents = function(x, highlight=c(), label="", tree_score=1) {
  if (purrr::is_empty(highlight)) highlight = get_unique_labels(x)

  # create an empty dataset with colnames
  edges = setNames(data.frame(matrix(ncol=3, nrow=0)), c("Parent", "Identity", "Label")) %>%
    tibble::as_tibble() %>%
    mutate(Parent=as.character(Parent), Identity=as.character(Identity), Label=as.character(Label))

  trees = get_trees(x, label)
  for (cluster in highlight) {
    tree = trees[[cluster]]

    if (!purrr::is_empty(tree))
      edges = edges %>%
        dplyr::add_row(
          tree[[tree_score]] %>%
            get_adj() %>%
            as.data.frame() %>%
            rownames_to_column(var="Label") %>%
            reshape2::melt(id="Label", variable.name="Identity") %>%
            filter(value==1) %>%
            filter(!Label %in% c("GL"), !Identity %in% c("GL", "P")) %>%
            dplyr::select(-value) %>%
            dplyr::mutate(Label=paste(cluster, Label, sep="."),
                          Identity=paste(cluster, Identity, sep="."),
                          Parent=cluster) %>%
            dplyr::mutate(Label=ifelse(grepl("P",Label), cluster, Label)) %>%
            tibble::as_tibble()
        )
  }
  return(edges)
}


# means format must be a dataframe with columns: labels, timepoints, lineage, mean_cov
get_muller_pop = function(x, map_tp_time=list("init"=0,"early"=60,"mid"=140,"late"=280),
                          mutations=FALSE, label="") {
  means = x %>% get_mean_long()

  if (mutations)
    # the means dataframe must contain also the subclones
    means = x %>%
      get_vaf_dataframe(label=label) %>%
      dplyr::select(theta, dplyr::contains("labels"), timepoints, lineage) %>%
      dplyr::select(-labels_init, -labels_viber) %>%
      tidyr::drop_na() %>%
      unique() %>%
      inner_join(means, by=c("labels", "timepoints", "lineage")) %>%
      mutate(mean_cov=theta*mean_cov) %>%
      dplyr::rename(parent=labels, labels=labels_mut) %>%
      dplyr::add_row( means %>% mutate(parent="P") )

  pop_df = means %>%
    format_means_df() %>%  # to format the dataframe with correct colnames ecc
    add_parent(x=x) %>%  # add common parent "P" data
    add_time_0(x=x, value="init") %>%
    convert_tp(mapping=map_tp_time) %>%  # convert timepoints to numeric values
    add_exp_fit_coeff(x=x) %>%
    dplyr::select(Identity, Generation, Lineage, Population, Frequency, dplyr::starts_with("lm"))

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


add_time_0 = function(pop_df, x=x, force=F, value="init") {
  n_tp = x %>% get_timepoints() %>% length()
  if (n_tp > 1 && !force) return(pop_df)

  print(value)

  ids = pop_df %>% filter(Identity!="P") %>% dplyr::pull(Identity) %>% unique()
  n_ids = length(ids)
  n_lins = x %>% get_lineages() %>% length()
  return(
    pop_df %>%
      dplyr::add_row(
        Identity=rep( ids, times = n_lins ),
        Population=rep( 0, times = n_ids * n_lins ),
        Frequency=rep( 0, times = n_ids * n_lins ),
        Generation=rep( value, times = n_ids * n_lins ),
        Lineage=rep( x$lineages, each = n_ids )
      ) %>%
      dplyr::add_row(
        Identity=rep( "P", times = n_lins ),
        Population=rep( 1, times = n_lins ),
        Frequency=rep( 1, times = n_lins ),
        Generation=rep( value, times = n_lins ),
        Lineage=x$lineages
      )
  )
}


convert_tp = function(pop_df, mapping=list("init"="0","early"="60","mid"="140","late"="280")) {
  if (is.null(mapping)) return(pop_df)

  mapping$init = "0"
  return(
    pop_df %>%
      mutate(Generation=mapping[Generation]) %>%
      mutate(Generation=unlist(Generation)) %>%
      mutate(Generation=as.numeric(Generation))
  )
}


add_exp_fit_coeff = function(pop_df, x) {
  if (pop_df$Generation %>% unique() %>% length() == 1) return(pop_df)
  return(
    pop_df %>%
      add_time_0(x=x, force=T, value=as.numeric(0)) %>%
      dplyr::mutate(Generation=ifelse(Generation=="init", 0, Generation)) %>%
      dplyr::group_by(Identity, Lineage) %>%
      dplyr::arrange(Generation, .by_group=T) %>%
      dplyr::mutate(lm_a=coef(lm(log1p(Population)~Generation))[1],
                    lm_r=coef(lm(log1p(Population)~Generation))[2]) %>%
      dplyr::ungroup() %>%
      dplyr::filter(Generation != 0)
  )
}


filter_muller_df = function(df, highlight=highlight) {
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


