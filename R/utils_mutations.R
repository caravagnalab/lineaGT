## Function to add the binomial theta, from the VIBER object, to the VAF dataframe
add_theta_to_vaf = function(vaf.df, x.muts.k, cluster) {
  theta = get_binomial_theta_cluster(x.muts.k, cluster)
  if (purrr::is_empty(theta))  # if there is no object fitted
    return(
      vaf.df %>%
        # mutate(theta_binom=as.numeric(vaf)/100, labels_mut=paste(labels, labels_binom, sep=".")) %>%
        mutate(theta_binom=as.numeric(vaf), labels_mut=paste(labels, labels_binom, sep=".")) %>%
        dplyr::select(-"labels_binom")
    )

  return(
    vaf.df %>%
      dplyr::mutate(labels_mut=paste(labels, labels_binom, sep=".")) %>%
      dplyr::inner_join(theta, by=c("labels_mut","labels","timepoints","lineage")) %>%
      dplyr::select(-"labels_binom")
  )
}


## Function to get and reshape the binomial theta from the VIBER object
get_binomial_theta_cluster = function(x.muts.k, cluster) {
  if (purrr::is_empty(x.muts.k))  # if no subclones, return an empty dataframe
    return(list())

  return(
    x.muts.k$theta_k %>%
      t() %>%
      as.data.frame() %>%
      dplyr::rename_with( ~ paste0("v.", .x)) %>%
      tibble::rownames_to_column(var="v_cluster") %>%
      mutate(labels_mut=paste(cluster, v_cluster, sep="."), v_cluster=NULL, labels=cluster) %>%
      tidyr::pivot_longer(cols=dplyr::starts_with("v"), names_to="v.tp.lin", values_to="theta_binom") %>%
      tidyr::separate("v.tp.lin", into=c("else","timepoints","lineage")) %>%
      mutate_tp(fn=as.integer) %>%
      dplyr::select(-"else")
  )
}


## Function to check for depth in mutations lower than "thr", substitute it with the mean coverage of the cluster,
## recompute "alt" and "ref" according to "vaf" values
check_dp = function(x, thr=10) {
  vaf.df = x %>% get_vaf_dataframe()
  means = x %>% get_mean_long()

  if (!"vaf" %in% colnames(vaf.df))
    vaf.df = vaf.df %>%
      dplyr::mutate(vaf=alt/dp)

  joined = dplyr::inner_join(vaf.df, means, by=c("labels", "timepoints", "lineage")) %>%

    # check 1 -> if "dp" is < "thr", set it as the mean coverage of the cluster
    dplyr::mutate(true_dp=dp) %>%
    dplyr::mutate(mean_cov=ceiling(mean_cov)) %>%
    dplyr::rowwise() %>% dplyr::mutate(dp=max(true_dp, mean_cov)) %>%
    # dplyr::group_by(labels) %>%
    # dplyr::rowwise() %>% dplyr::mutate(dp=replace(dp, any(dp < thr) & (dp < thr), mean_cov)) %>%
    # dplyr::ungroup() %>%

    # check 2 -> if the dp still is < "thr", set it to the mean depth of the other mutations
    dplyr::group_by(labels) %>%
    dplyr::rowwise() %>% dplyr::mutate(dp=replace(dp, dp < thr, mean(dp) %>% as.integer())) %>%
    dplyr::ungroup() %>%

    # last check -> if still < "thr" set it to thr
    dplyr::mutate(dp=replace(dp, dp < thr, thr)) %>%

    # correct the variant allele reads
    # dplyr::mutate(alt=ceiling(vaf/100*dp)) %>%
    dplyr::mutate(alt=ceiling(vaf*dp)) %>%

    dplyr::select(dp, alt, vaf, mutation, IS, lineage, timepoints, true_dp, dplyr::starts_with("ref"),
                  dplyr::starts_with("labels"), dplyr::starts_with("theta"), dplyr::starts_with("pi"))

    # # check 1 -> if "dp" is < "thr", set it as the mean coverage of the cluster
    # dplyr::mutate(true_dp=dp) %>%
    # dplyr::mutate(mean_cov=ceiling(mean_cov)) %>%
    # dplyr::rowwise() %>% dplyr::mutate(dp=replace(dp, dp < thr, mean_cov)) %>%  # replace values < thr to mean_cov
    #
    # # check 2 -> if the dp still is < "thr", set it to the mean depth of the other mutations
    # dplyr::group_by(labels) %>%
    # dplyr::mutate(dp=replace(dp, dp < thr, mean(dp) %>% as.integer())) %>%
    # dplyr::ungroup() %>%
    #
    # # last check -> if still < "thr" set it to thr
    # dplyr::mutate(dp=replace(dp, dp < thr, thr)) %>%
    #
    # # correct the variant allele reads
    # dplyr::mutate(alt=ceiling(vaf/100*dp)) %>%
    #
    # dplyr::select(dp, alt, vaf, mutation, IS, lineage, timepoints, true_dp, dplyr::starts_with("ref"),
    #               dplyr::starts_with("labels"), dplyr::starts_with("theta"), dplyr::starts_with("pi"))

  return(add_vaf(x, joined))
}


# Function to check if dimensions found in the coverage dataset (as combinations
# of lineages and timepoints) are missing in the mutation data
check_vaf_dimensions = function(vaf.df, x) {
  vaf.df = vaf.df %>%
    dplyr::ungroup() %>%
    long_to_wide_muts() %>%
    wide_to_long_muts()

  vaf.dims = vaf.df %>%
    dplyr::group_by(lineage, timepoints) %>%
    dplyr::summarise(nn=dplyr::n(), .groups="keep") %>%
    dplyr::mutate(dimensions=paste("cov",timepoints,lineage,sep=".")) %>%
    dplyr::pull(dimensions)

  cov.dims = x %>% get_dimensions()

  missing.vals = setdiff(cov.dims, vaf.dims) %>%
    reshape2::melt() %>%
    separate(value, into=c("else","timepoints","lineage")) %>%
    dplyr::mutate("else"=NULL)

  if (nrow(missing.vals) == 0) return( vaf.df )

  # try(expr = {
  #   missing.vals = missing.vals %>%
  #     dplyr::mutate(timepoints=as.integer(timepoints))
  # })

  return(
    vaf.df %>%
      dplyr::add_row(
        missing.vals %>%
          mutate_tp(fn=as.integer, colnm="timepoints") %>%
          dplyr::mutate(mutation=vaf.df[1,] %>% dplyr::pull(mutation),
                        IS=vaf.df[1,] %>% dplyr::pull(IS))
        ) %>%
      replace(is.na(.), 0) %>%
      long_to_wide_muts() %>%
      wide_to_long_muts()
    )
}


# returns the object x with the annotated vaf dataframe
annotate_vaf_df = function(x, vaf.df, min_frac=0) {
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
                  dplyr::starts_with("labels"),
                  dplyr::starts_with("theta"),
                  dplyr::starts_with("pi"),
                  mutation, IS, lineage, timepoints)

  vaf.ann = dplyr::inner_join(vaf.df_filt, dataframe, by=c("IS", "lineage", "timepoints"))

  return(add_vaf(x, vaf.ann))
}


get_input_viber = function(x) {

  vaf.df_wide = x %>%
    get_vaf_dataframe() %>%
    long_to_wide_muts()

  trials = vaf.df_wide %>%
    dplyr::select(starts_with("dp."), labels) %>%
    rename_with(.fn = ~str_replace_all(.x, "dp.", ""))

  successes = vaf.df_wide %>%
    dplyr::select(starts_with("alt."), labels) %>%
    rename_with(.fn = ~str_replace_all(.x, "alt.", ""))

  return(list("successes"=successes, "trials"=trials, "vaf.df"=vaf.df_wide))

  # alpha_0 = x %>%
  #   get_vaf_dataframe() %>%
  #   group_by(labels, lineage, timepoints) %>%
  #   dplyr::summarise(mean_vaf=mean(vaf)) %>%
  #   dplyr::mutate(alpha_0=ifelse(mean_vaf==0, 0.01,1)) %>%
  #   dplyr::select(-mean_vaf) %>%
  #   tidyr::pivot_wider(names_from=c("timepoints","lineage"), values_from="alpha_0", names_sep=".")

  # return(list("successes"=successes, "trials"=trials, "vaf.df"=vaf.df_wide, "alpha_0"=alpha_0))
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
           dplyr::filter(labels==cluster))

         # "alpha_0"=input$alpha_0 %>%
         #   dplyr::filter(labels==cluster))
  )
}


get_data_annotation = function(k) {
  data_annotations = data.frame(gene=paste0("G", 1:k), driver=FALSE)
  data_annotations$driver[sample(1:nrow(data_annotations), 1)] = TRUE
  return(data_annotations)
}


