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


get_mean_wide = function(x, arrange=T) {
  return(
    x %>%
      get_mean_long() %>%
      group_by(labels) %>%
      dplyr::arrange(desc(mean_cov)) %>%
      pivot_wider(id_cols="labels", names_from=c("timepoints","lineage"), values_from="mean_cov")
  )
}


check_dp = function(x, thr=5, label="") {
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
    dplyr::mutate(dp=ifelse(dp < thr, mean(dp) %>% as.integer(), dp)) %>%
    dplyr::ungroup() %>%

    dplyr::select(-original_dp)

  return(add_vaf(x, joined, label))
}


# returns the object x with the annotated vaf dataframe
annotate_vaf_df = function(x, vaf.df, min_frac=0, label="") {
  highlight = get_highlight(x, min_frac=min_frac)

  dataframe = x %>%
    get_cov_dataframe() %>%
    dplyr::filter(labels %in% highlight)
  IS_keep = dataframe$IS %>% unique()

  vaf.df_filt = vaf.df %>%
    dplyr::filter(IS %in% IS_keep,
                  lineage %in% (x %>% get_lineages()),
                  timepoints %in% (x %>% get_timepoints())) %>%
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
  )
}


get_data_annotation = function(k) {
  data_annotations = data.frame(gene=paste0("G", 1:k), driver=FALSE)
  data_annotations$driver[sample(1:nrow(data_annotations), 1)] = TRUE
  return(data_annotations)
}


# # As input a mvnmm object with already a viber_run performed
# get_binomial_theta = function(x, label="") {
#   x.muts = x %>% get_viber_run(label=label)
#   theta = data.frame()
#
#   for(cluster in x.muts %>% names()) {
#     if (!purrr::is_empty(x.muts[[cluster]])) {
#       theta.k = get_binomial_theta_cluster(x.muts[[cluster]], cluster)
#       theta = rbind(theta, theta.k)
#     }
#   }
#
#   if (!purrr::is_empty(theta))
#     return(
#       theta %>%
#         tidyr::pivot_longer(cols=starts_with("v."), names_to="v.timepoints.lineage", values_to="theta") %>%
#         tidyr::separate(v.timepoints.lineage, into=c("else","timepoints","lineage")) %>%
#         dplyr::mutate("else"=NULL, theta=theta*100)
#     )
# }


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
