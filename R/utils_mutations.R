# As input a mvnmm object with already a viber_run performed
get_binomial_theta = function(x) {
  viber_fits = x$viber_run
  theta = data.frame()
  for(cluster in viber_fits %>% names()) {
    if (!purrr::is_empty(viber_fits[[cluster]])) {
      df_k = viber_fits[[cluster]]$theta_k %>% t %>% as.data.frame() %>%
        dplyr::rename_with( ~ paste0("v.", .x)) %>% tibble::rownames_to_column(var="v_cluster") %>%
        mutate(labels_mut=paste(cluster, v_cluster, sep="."), v_cluster=NULL, labels=cluster)
      theta = rbind(theta, df_k)
    }
  }
  return(
    theta %>%
      tidyr::pivot_longer(cols=starts_with("v."), names_to="v.timepoints.lineage", values_to="theta") %>%
      separate(v.timepoints.lineage, into=c("else","timepoints","lineage")) %>%
      mutate("else"=NULL)
  )
}



# Functions used to obtain and reshape some datasets
# works with our vaf files
vaf_df_from_file = function(vaf_file) {
  vaf_df = read.csv(vaf_file)
  try(expr = { vaf_df = vaf_df  %>%
    dplyr::rename_with(.cols=all_of(dplyr::starts_with("dp_")),
                       .fn=~paste0("cov_", str_replace_all(.x,"dp_",""))) }, silent = T)

  vaf_df = vaf_df %>%
    tidyr::pivot_longer(cols=starts_with("dp.ref.alt"), names_to="timepoint", values_to="dp:ref:alt") %>%
    mutate(timepoint=stringr::str_replace_all(timepoint, "dp.ref.alt_", "")) %>%
    filter(!timepoint %in% c("over","steady")) %>%
    separate("dp:ref:alt", into=c("dp", "ref", "alt"), sep=":") %>%
    mutate(ref=as.integer(ref), alt=as.integer(alt), dp=ref+alt) %>%
    tidyr::pivot_wider(values_from=c("dp","ref","alt"), names_from="timepoint", values_fn=as.integer)

  vaf_df = vaf_df %>%
    dplyr::select(-starts_with("cov")) %>%
    tidyr::pivot_longer(cols=c(starts_with("alt"),starts_with("dp"),starts_with("vaf"),starts_with("ref"))) %>%
    separate(name, into=c("type","timepoints")) %>%
    tidyr::pivot_wider(names_from = "type", values_from = "value")

  return(vaf_df)
}


check_dp = function(vaf.df, x) {
  # vaf_df is the dataframe from vaf_df_from_file
  # returns a dataframe with same structure but dp_early/mid/late equals the
  # mean of the other timepoints when 0
  # mean coverage for the corresponding cluster when 0 again
  return(
    vaf.df %>%
      group_by(lineage, mutation, IS) %>%
      mutate(dp=ifelse(dp==0, ceiling(mean(dp)), dp)) %>%
      ungroup() %>%
      mutate(dp=ifelse(dp==0, ceiling(get_mean(x)[labels,paste("cov",timepoints,lineage,sep=".")]), dp))
  )
}


annotate_vaf_df = function(vaf.df, x, min_frac=0) {
  ll = x$lineages
  tp = x$timepoints
  highlight = select_relevant_clusters(x, min_frac=min_frac)
  dataframe = x %>% get_dataframe() %>% filter(labels %in% highlight)
  IS_keep = dataframe$IS %>% unique()

  vaf.df_filt = vaf.df %>%
    filter(IS %in% IS_keep, lineage %in% ll, timepoints %in% tp) %>%
    dplyr::select(alt, dp, vaf, mutation, IS, lineage, timepoints)

  return(dplyr::inner_join(vaf.df_filt, dataframe, by=c("IS", "lineage", "timepoints")))
}


# Function to get from a vaf dataframe obtained by get_vaf_df() the input for a VIBER run
get_input_viber = function(vaf.df, x) {
  vaf.df_wide = check_dp(vaf.df, x) %>% long_to_wide_muts()

  trials = vaf.df_wide %>% dplyr::select(starts_with("dp"), labels) %>%
    rename_with(.fn=~str_replace_all(.x,"dp.",""))

  successes = vaf.df_wide %>% dplyr::select(starts_with("alt"), labels) %>%
    rename_with(.fn=~str_replace_all(.x,"alt.",""))

  return(list("successes"=successes, "trials"=trials, "vaf.df"=vaf.df_wide))
}



wide_to_long_muts = function(vaf.df) {
  return(
    vaf.df %>%
      tidyr::pivot_longer(cols=c(starts_with("alt"), starts_with("dp"), starts_with("vaf")),
                          names_to="type.timepoints.lineage") %>%
      separate(type.timepoints.lineage, into=c("type", "timepoints", "lineage")) %>%
      tidyr::pivot_wider(names_from="type", values_from="value")
  )
}



long_to_wide_muts = function(vaf.df) {
  return(
    vaf.df %>%
      dplyr::select(alt, dp, vaf, timepoints, lineage, IS, mutation, starts_with("labels")) %>%
      tidyr::pivot_wider(names_from=c("timepoints","lineage"), names_sep=".",
                         values_from=c("alt","dp","vaf"))
  )
}
