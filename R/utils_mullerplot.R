format_means_df = function(mean_df) {
  return(
    mean_df %>%
      dplyr::rename(Identity=labels, Generation=timepoints, Lineage=lineage, Population=pop) %>%
      dplyr::mutate(Identity=as.character(Identity)) %>%
      dplyr::mutate(Population=ifelse(Population==0, 0.001, Population)) %>%
      dplyr::arrange(Identity, Lineage, Generation)
  )
}


add_parent = function(pop_df, x) {
  return(
    pop_df %>%
      dplyr::add_row(
        Identity=rep( "P", times = x$data.shape[2] ),
        Population=rep( 0, times = x$data.shape[2] ),
        Pop.plot=rep( 0, times = x$data.shape[2] ),
        Frequency=rep( 0, times = x$data.shape[2] ),
        theta_binom=rep( 0, times = x$data.shape[2] ),
        Generation=rep( get_timepoints(x), times = get_lineages(x) %>% length() ),
        Lineage=rep( x$lineages, each = get_timepoints(x) %>% length() )
        ) %>%
      dplyr::arrange(Identity, Lineage, Generation)
    )
}


add_time_0 = function(pop_df,
                      x,
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
        Pop.plot=rep( 0, times = n_ids * n_lins ),
        Frequency=rep( 0, times = n_ids * n_lins ),
        theta_binom=rep( 0, times = n_ids * n_lins ),
        Generation=rep( value, times = n_ids * n_lins ),
        Lineage=rep( x$lineages, each = n_ids )
      ) %>%
      dplyr::add_row(
        Identity=rep( "P", times = n_lins ),
        Population=rep( 1, times = n_lins ),
        Pop.plot=rep( 1, times = n_lins ),
        Frequency=rep( 1, times = n_lins ),
        theta_binom=rep( 1, times = n_lins ),
        Generation=rep( value, times = n_lins ),
        Lineage=x$lineages
      ) %>%
      dplyr::arrange(Identity, Lineage, Generation)
  )
}


convert_tp = function(pop_df,
                      timepoints_to_int) {
  # if (is.null(timepoints_to_int)) return(pop_df)
  if(is.numeric(pop_df$Generation)) pop_df = pop_df %>% dplyr::mutate(Generation=as.character(Generation))
  return(
    pop_df %>%
      dplyr::mutate(Generation=timepoints_to_int[Generation])
  )
}


filter_muller_df = function(pop_df, highlight=highlight) {
  return(
    pop_df %>%
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
    dplyr::ungroup()
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
    dplyr::mutate(Frequency=NA,
           Group_id="___special_emptya",
           Unique_id=paste0("___special_emptya_",Generation))

  mullerdf = new_rows1 %>%
    dplyr::bind_rows(mullerdf) %>%
    dplyr::arrange(Generation) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(new_rows2) %>%
    dplyr::arrange(Generation) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Generation, Lineage) %>%
    dplyr::mutate(Frequency=Population/sum(Population)) %>%
    dplyr::ungroup()

  mullerdf$Group_id = factor(mullerdf$Group_id,
                             levels=rev(unlist(
                               as.data.frame(mullerdf %>%
                                               dplyr::filter_(~Generation == max(Generation)) %>%
                                               dplyr::select_(~Group_id)),
                               use.names=FALSE) %>% unique()))
  return(mullerdf)
}


estimate_mean_ISs = function(x) {
  cov.df = x %>% get_cov_dataframe() %>% long_to_wide_cov() %>% dplyr::select(dplyr::starts_with("cov"))

  cums = apply(cov.df, 2, function(x) return(ecdf(x)))
  qntls = apply(cov.df, 2, function(x)
    return(max(1, quantile(x, probs=0.95) %>% setNames(NULL)))) %>% unlist()

  qntls.df = data.frame(qntls) %>%
    tibble::rownames_to_column() %>%
    tidyr::separate(rowname, into=c("else","timepoints","lineage")) %>%
    dplyr::select(-"else") %>%
    tibble::as_tibble() %>%
    mutate_tp(fn=as.integer)

  mean.long = x %>% get_mean_long()

  cov.qntls = dplyr::inner_join(mean.long, qntls.df, by=c("timepoints","lineage"))

  keep = cov.qntls %>%
    dplyr::group_by(labels) %>%
    dplyr::filter(any(mean_cov>qntls)) %>%
    dplyr::pull(labels) %>% unique()

  if (length(keep) == 0)
    return(get_ISs(x) %>% mean() %>% ceiling())

  # fixx = setdiff(get_unique_labels(x), keep)
  mean_ISs = get_ISs(x)[keep] %>% mean() %>% ceiling()

  return(mean_ISs)
}


estimate_n_pops = function(x, highlight=c()) {
  n_ISs = get_ISs(x)
  mean_ISs = estimate_mean_ISs(x)

  n_pops = sapply(get_unique_labels(x), function(cls) return(max(1, round(n_ISs[[cls]] / mean_ISs))) )

  return(n_pops)
}


