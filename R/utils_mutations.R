# Functions used to obtain and reshape some datasets
#
# NOTE: get_vaf_df returns a dataset with mutation,IS,lineage,timepoints,dp,ref,alt,vaf,... columns
# usage:
# vaf_df = vaf_df_from_file(vaf_file)
# vaf_df = annotate_vaf_df(vaf_df, obj, min_ccf=0.07)
# obj = add_vaf(obj, vaf_df)


vaf_df_from_file = function(vaf_file) {
  ## input is the file with the vaf, having a column per timepoint named dp.ref.alt
  ## final dataframe will have:
  ## mutation, lineage, IS, mut_type
  ## vaf_early/mid/late,
  ## cov_early/mid/late,
  ## dp_early/mid/late,
  ## ref_early/mid/late
  ## alt_early/mid/late
  vaf_df = read.csv(vaf_file)
  try(expr = { vaf_df = vaf_df  %>%
    dplyr::rename_with(.cols=all_of(dplyr::starts_with("dp_")),
                       .fn=~paste0("cov_", str_replace_all(.x,"dp_",""))) }, silent = T)

  vaf_df = vaf_df %>%
    tidyr::pivot_longer(cols=starts_with("dp.ref.alt"), names_to="timepoint", values_to="dp:ref:alt") %>%
    mutate(timepoint=stringr::str_replace_all(timepoint, "dp.ref.alt_", "")) %>%
    separate("dp:ref:alt", into=c("dp", "ref", "alt"), sep=":") %>%
    tidyr::pivot_wider(values_from=c("dp","ref","alt"), names_from="timepoint", values_fn=as.integer) %>%
    update_trials()

  return(vaf_df)
}


update_trials = function(vaf_df) {
  # vaf_df is the dataframe from vaf_df_from_file
  # returns a dataframe with same structure but dp_early/mid/late equals the
  # mean of the other timepoints when 0

  vaf_df_new = vaf_df %>%
    tidyr::pivot_longer(cols=c(starts_with("dp")), names_to=c("timepoints"), values_to="dp",
                        values_transform=list(dp=as.integer)) %>%
    group_by(lineage, mutation, IS) %>%
    mutate(dp=ifelse(dp==0, mean(dp), dp)) %>% dplyr::ungroup() %>%
    tidyr::pivot_wider(values_from="dp", names_from="timepoints")

  return(vaf_df_new)
}


annotate_vaf_df = function(vaf_df, obj, min_ccf=0) {
  lineages = obj$lineages
  highlight = select_relevant_clusters(obj, min_ccf=min_ccf)

  dataframe = obj$dataframe %>% filter(labels %in% highlight) %>%
    tidyr::pivot_longer(cols=starts_with("cov"), names_to="cov_timepoints_lineage",
                        values_to="cov", values_transform=list(cov=as.integer)) %>%
    separate(cov_timepoints_lineage, into=c("cc", "timepoints", "lineage")) %>%
    tidyr::pivot_wider(names_from=c("cc","timepoints"), values_from="cov")

  IS_keep = dataframe$IS

  vaf_df_filt = vaf_df %>% filter(IS %in% IS_keep, lineage %in% lineages) %>%
    dplyr::select(starts_with("alt"), starts_with("dp"), starts_with("vaf"), mutation, IS, lineage)

  vaf_annotated = dplyr::inner_join(vaf_df_filt, dataframe, by=c("IS", "lineage"))

  return(vaf_annotated)
}


# Function to get from a vaf dataframe obtained by get_vaf_df() the input for a VIBER run
get_input_viber = function(vaf_df) {
  vaf_df_wide = vaf_df %>%
    tidyr::pivot_wider(values_from=c(starts_with("dp"), starts_with("alt"), starts_with("ref"),
                                     starts_with("vaf"), starts_with("cov")),
                       names_from=lineage, values_fn=as.numeric)
  trials = vaf_df_wide %>% dplyr::select(starts_with("dp"), labels) %>%
    rename_with(.fn=~str_replace_all(.x,"dp_",""))
  successes = vaf_df_wide %>% dplyr::select(starts_with("alt"), labels) %>%
    rename_with(.fn=~str_replace_all(.x,"alt_",""))

  return(list("successes"=successes, "trials"=trials, "vaf_df"=vaf_df_wide))
}


# As input a mvnmm object with already the fitted viber df
vaf_dataframe = function(obj) {
  joined = obj$dataframe_vaf %>% tidyr::unite(col="labels_mut", c("labels", "labels_viber"), sep=".", remove=F)
  dp = joined %>%
    reshape2::melt(id=c("IS", "labels", "mutation", "experiment", "pi_viber", "labels_mut", "labels_viber")) %>%
    tidyr::separate(variable, into=c("time_lin", "dp_alt"), sep="[.]") %>%
    filter(dp_alt == "dp") %>% mutate(dp=value, value=NULL) %>% dplyr::select(-dp_alt)
  alt = joined %>%
    reshape2::melt(id=c("IS", "labels", "mutation", "experiment", "pi_viber", "labels_mut", "labels_viber")) %>%
    tidyr::separate(variable, into=c("time_lin", "dp_alt"), sep="[.]") %>%
    filter(dp_alt == "alt") %>% mutate(alt=value, value=NULL) %>% dplyr::select(-dp_alt)
  vaf = dplyr::inner_join(dp, alt) %>% mutate(vaf=alt/dp*100) %>% mutate(vaf=ifelse(is.na(vaf), 0, vaf)) %>%
    tidyr::separate("time_lin", into=c("timepoint", "lineage")) %>% dplyr::select(-"dp", -"alt") %>%
    tidyr::pivot_wider(names_from=c("timepoint"), values_from="vaf", names_sep=".", names_prefix="vaf.")

  return(vaf)
}


# As input a mvnmm object with already a viber_run performed
get_binomial_theta = function(obj) {
  viber_fits = obj$viber_run
  theta = data.frame()
  for(cluster in viber_fits %>% names()) {
    if (class(viber_fits[[cluster]]) != "character") {
      df_k = viber_fits[[cluster]]$theta_k %>% t %>% as.data.frame() %>%
        dplyr::rename_with( ~ paste0("vaf.", .x)) %>% tibble::rownames_to_column(var="v_cluster") %>%
        mutate(labels_mut=paste(cluster, v_cluster, sep="."), v_cluster=NULL, labels=cluster)
      theta = rbind(theta, df_k)
    }
  }
  theta = theta %>% reshape2::melt() %>% mutate(value=value*100) %>%
    tidyr::separate(variable, into=c("timepoint", "lineage"), sep="_") %>%
    tidyr::pivot_wider(names_from="timepoint", values_from="value")
  return(theta)
}


