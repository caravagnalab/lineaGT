# As input a mvnmm object with already a viber_run performed
get_binomial_theta = function(x, label="") {
  viber_fits = x %>% get_viber_run(label=label)
  theta = data.frame()
  for(cluster in viber_fits %>% names()) {
    if (!purrr::is_empty(viber_fits[[cluster]])) {
      df_k = get_binomial_theta_cluster(viber_fits[[cluster]], cluster)
      theta = rbind(theta, df_k)
    }
  }
  return(
    theta %>%
      tidyr::pivot_longer(cols=starts_with("v."), names_to="v.timepoints.lineage", values_to="theta") %>%
      separate(v.timepoints.lineage, into=c("else","timepoints","lineage")) %>%
      mutate("else"=NULL, theta=theta*100)
  )
}


get_binomial_theta_cluster = function(viber_fit, cluster) {
  return(
    viber_fit$theta_k %>%
      t() %>%
      as.data.frame() %>%
      dplyr::rename_with( ~ paste0("v.", .x)) %>%
      tibble::rownames_to_column(var="v_cluster") %>%
      mutate(labels_mut=paste(cluster, v_cluster, sep="."), v_cluster=NULL, labels=cluster)
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


check_dp = function(x, thr=5, label="") {
  vaf.df = x %>% get_vaf_dataframe(label)
  means = x %>% get_mean_long()

  joined = dplyr::inner_join(vaf.df, means, by=c("labels", "timepoints", "lineage")) %>%
    dplyr::mutate(original_dp=dp) %>%
    dplyr::mutate(dp=ceiling(mean_cov)) %>%
    dplyr::rowwise() %>% dplyr::mutate(dp=max(dp,original_dp)) %>%
    dplyr::mutate(alt=ceiling(vaf/100*dp)) %>%
    dplyr::group_by(labels) %>%
    dplyr::mutate(dp=ifelse(dp < thr, mean(dp) %>% as.integer(), dp)) %>%
    dplyr::ungroup()

  return(add_vaf(x, joined, label))
}


# returns the object x with the annotated vaf dataframe
annotate_vaf_df = function(x, vaf.df, min_frac=0, label="") {
  ll = x %>% get_lineages()
  tp = x %>% get_timepoints()
  highlight = select_relevant_clusters(x, min_frac=min_frac)

  dataframe = x %>%
    get_cov_dataframe() %>%
    dplyr::filter(labels %in% highlight)
  IS_keep = dataframe$IS %>% unique()

  vaf.df_filt = vaf.df %>%
    dplyr::filter(IS %in% IS_keep, lineage %in% ll, timepoints %in% tp) %>%
    dplyr::select(alt, dp, vaf, mutation, IS, lineage, timepoints)

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
    dplyr::select(starts_with("dp"), labels) %>%
    rename_with(.fn = ~str_replace_all(.x, "dp.", ""))

  successes = vaf.df_wide %>%
    dplyr::select(starts_with("alt"), labels) %>%
    rename_with(.fn = ~str_replace_all(.x, "alt.", ""))
  return(list("successes"=successes, "trials"=trials, "vaf.df"=vaf.df_wide))
}


# filter_muts = function(vaf.df) {
#   vaf.df %>%
#     group_by(mutation, lineage, labels) %>%
#     dplyr::mutate(vaf_diff=vaf-dplyr::lag(vaf), vaf_diff=ifelse(is.na(vaf_diff),0,vaf_diff)) %>%
#     dplyr::filter(any(abs(vaf_diff)>40))
# }


# # Functions used to obtain and reshape some datasets
# # works with our vaf files
# vaf_df_from_file = function(vaf_file) {
#   vaf_df = read.csv(vaf_file)
#   try(expr = { vaf_df = vaf_df  %>%
#     dplyr::rename_with(.cols=all_of(dplyr::starts_with("dp_")),
#                        .fn=~paste0("cov_", str_replace_all(.x,"dp_",""))) }, silent = T)
#
#   vaf_df = vaf_df %>%
#     tidyr::pivot_longer(cols=starts_with("dp.ref.alt"), names_to="timepoint", values_to="dp:ref:alt") %>%
#     mutate(timepoint=stringr::str_replace_all(timepoint, "dp.ref.alt_", "")) %>%
#     filter(!timepoint %in% c("over","steady")) %>%
#     separate("dp:ref:alt", into=c("dp", "ref", "alt"), sep=":") %>%
#     mutate(ref=as.integer(ref), alt=as.integer(alt), dp=ref+alt) %>%
#     tidyr::pivot_wider(values_from=c("dp","ref","alt"), names_from="timepoint", values_fn=as.integer)
#
#   vaf_df = vaf_df %>%
#     dplyr::select(-starts_with("cov")) %>%
#     tidyr::pivot_longer(cols=c(starts_with("alt"),starts_with("dp"),starts_with("vaf"),starts_with("ref"))) %>%
#     separate(name, into=c("type","timepoints")) %>%
#     tidyr::pivot_wider(names_from = "type", values_from = "value")
#
#   return(vaf_df)
# }
