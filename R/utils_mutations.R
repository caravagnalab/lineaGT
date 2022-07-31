add_theta_to_vaf = function(vaf.df, x.muts.k, cluster, label="") {
  theta = get_binomial_theta_cluster(x.muts.k, cluster)
  if (purrr::is_empty(theta))
    return(
      vaf.df %>%
        mutate(theta=NA, labels_mut=paste(labels, labels_viber, sep="."))
    )
  return(
    vaf.df %>%
      dplyr::mutate(labels_mut=paste(labels, labels_viber, sep=".")) %>%
      dplyr::inner_join(theta, by=c("labels_mut","labels","timepoints","lineage"))
  )
}


get_binomial_theta_cluster = function(x.muts.k, cluster) {
  if (purrr::is_empty(x.muts.k))
    # if no subclones, return an empty dataframe
    return(list())
  return(
    x.muts.k$theta_k %>%
      t() %>%
      as.data.frame() %>%
      dplyr::rename_with( ~ paste0("v.", .x)) %>%
      tibble::rownames_to_column(var="v_cluster") %>%
      mutate(labels_mut=paste(cluster, v_cluster, sep="."), v_cluster=NULL, labels=cluster) %>%
      tidyr::pivot_longer(cols=dplyr::starts_with("v"), names_to="v.tp.lin", values_to="theta") %>%
      tidyr::separate("v.tp.lin", into=c("else","timepoints","lineage")) %>%
      dplyr::select(-"else")
  )
}


get_mean_long = function(x) {
  return(
    x %>%
      get_mean() %>%
      as.data.frame() %>% tibble::rownames_to_column(var="labels") %>%
      dplyr::mutate(labels=factor(labels, levels=unique(labels))) %>%
      tidyr::pivot_longer(cols=starts_with("cov"), names_to="cov.timepoints.lineage", values_to="mean_cov") %>%
      tidyr::separate(cov.timepoints.lineage, into=c("else", "timepoints", "lineage"), sep="[.]") %>%
      dplyr::mutate("else"=NULL)
  )
}


get_variance_long = function(x) {
  return(
    x %>%
      get_sigma() %>%
      as.data.frame() %>% tibble::rownames_to_column(var="labels") %>%
      dplyr::mutate(labels=factor(labels, levels=unique(labels))) %>%
      tidyr::pivot_longer(cols=starts_with("cov"), names_to="cov.timepoints.lineage", values_to="sigma") %>%
      tidyr::separate(cov.timepoints.lineage, into=c("else", "timepoints", "lineage"), sep="[.]") %>%
      dplyr::mutate("else"=NULL)
  )
}


get_mean_wide = function(x, arrange=T) {
  return(
    x %>%
      get_mean_long() %>%
      group_by(labels) %>%
      dplyr::arrange(desc(mean_cov)) %>%
      pivot_wider(id_cols="labels", names_from=c("timepoints","lineage"), values_from="mean_cov")
  )
}


check_dp = function(x, thr=10, label="") {
  vaf.df = x %>% get_vaf_dataframe(label=label)
  means = x %>% get_mean_long()

  joined = dplyr::inner_join(vaf.df, means, by=c("labels", "timepoints", "lineage")) %>%

    # set the depth `dp` as the maximum value among the called depth and the mean coverage
    dplyr::mutate(original_dp=dp) %>%
    dplyr::mutate(dp=ceiling(mean_cov)) %>%
    dplyr::rowwise() %>% dplyr::mutate(dp=max(dp,original_dp)) %>%
    dplyr::mutate(alt=ceiling(vaf/100*dp)) %>%

    # check the minimum depth to be at least `thr`, otherwise
    dplyr::group_by(labels) %>%
    dplyr::mutate(dp=replace(dp, dp < thr, mean(dp) %>% as.integer())) %>%
    dplyr::ungroup() %>%

    dplyr::select(-original_dp, dp, alt, ref, vaf, mutation, IS, lineage, timepoints, labels)

  return(add_vaf(x, joined, label))
}


# Function to check if dimensions found in the coverage dataset (as combinations
# of lineages and timepoints) are missing in the mutation data
check_vaf_dimensions = function(vaf.df, x) {
  vaf.df = vaf.df %>%
    ungroup() %>%
    long_to_wide_muts() %>%
    wide_to_long_muts()

  vaf.dims = vaf.df %>%
    group_by(lineage, timepoints) %>%
    dplyr::summarise(nn=dplyr::n()) %>%
    mutate(dimensions=paste("cov",timepoints,lineage,sep=".")) %>%
    dplyr::pull(dimensions)
  cov.dims = x %>% get_dimensions()

  missing = setdiff(cov.dims, vaf.dims) %>%
    reshape2::melt() %>%
    separate(value, into=c("else","timepoints","lineage")) %>%
    dplyr::mutate("else"=NULL)

  if (purrr::is_empty(missing))
    return(
      vaf.df
    )

  return(
    vaf.df %>%
      dplyr::add_row(
        missing %>%
          dplyr::mutate(mutation=vaf.df[1,] %>% dplyr::pull(mutation),
                        IS=vaf.df[1,] %>% dplyr::pull(IS),
                        dp=0, dp_locus=0, alt=0, ref=0, vaf=0.0)
        ) %>%
      long_to_wide_muts() %>%
      wide_to_long_muts()
    )
}


# returns the object x with the annotated vaf dataframe
annotate_vaf_df = function(x, vaf.df, min_frac=0, label="") {
  highlight = get_highlight(x, min_frac=min_frac)

  vaf.df = vaf.df %>% check_vaf_dimensions(x=x)

  dataframe = x %>%
    get_cov_dataframe() %>%
    dplyr::filter(labels %in% highlight)
  IS_keep = dataframe$IS %>% unique()

  vaf.df_filt = vaf.df %>%
    dplyr::filter(IS %in% IS_keep,
                  lineage %in% (x %>% get_lineages()),
                  timepoints %in% (x %>% get_timepoints())) %>%
    dplyr::select(dplyr::starts_with("alt"),
                  dplyr::starts_with("ref"),
                  dplyr::starts_with("dp"),
                  dplyr::starts_with("vaf"),
                  mutation, IS, lineage, timepoints)

  vaf.ann = dplyr::inner_join(vaf.df_filt, dataframe, by=c("IS", "lineage", "timepoints"))

  return(add_vaf(x, vaf.ann, label))
}


# Function to get from a vaf dataframe obtained by vaf_df_from_file() the input for a VIBER run
get_input_viber = function(x, lineages=c(), label="") {
  if (purrr::is_empty(lineages)) lineages = x %>% get_lineages()

  vaf.df_wide = x %>%
    get_vaf_dataframe(label) %>%
    dplyr::filter(lineage%in%lineages) %>%
    long_to_wide_muts()

  trials = vaf.df_wide %>%
    dplyr::select(starts_with("dp."), labels) %>%
    rename_with(.fn = ~str_replace_all(.x, "dp.", ""))

  successes = vaf.df_wide %>%
    dplyr::select(starts_with("alt."), labels) %>%
    rename_with(.fn = ~str_replace_all(.x, "alt.", ""))

  alpha_0 = x %>%
    get_vaf_dataframe(label) %>%
    group_by(labels, lineage, timepoints) %>%
    dplyr::summarise(mean_vaf=mean(vaf)) %>%
    dplyr::mutate(alpha_0=ifelse(mean_vaf==0, 0.01,1)) %>%
    dplyr::select(-mean_vaf) %>%
    tidyr::pivot_wider(names_from=c("timepoints","lineage"), values_from="alpha_0", names_sep=".")

  return(list("successes"=successes, "trials"=trials, "vaf.df"=vaf.df_wide, "alpha_0"=alpha_0))
}


filter_viber_input = function(input, cluster) {
  return(
    list("successes"=input$successes %>%
           dplyr::filter(labels==cluster) %>%
           dplyr::select(-labels),

         "trials"=input$trials %>%
           dplyr::filter(labels==cluster) %>%
           dplyr::select(-labels),

         "vaf.df"=input$vaf.df %>%
           dplyr::filter(labels==cluster),

         "alpha_0"=input$alpha_0 %>%
           dplyr::filter(labels==cluster))
  )
}


get_data_annotation = function(k) {
  data_annotations = data.frame(gene=paste0("G", 1:k), driver=FALSE)
  data_annotations$driver[sample(1:nrow(data_annotations), 1)] = TRUE
  return(data_annotations)
}


