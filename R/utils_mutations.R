#' Functions used to obtain and reshape some datasets
#'
#' NOTE: get_vaf_df returns a dataset with mutation,IS,lineage,timepoints,dp,ref,alt,vaf,... columns
#'
#' @importFrom tidyr separate unite pivot_wider pivot_longer
#' @importFrom dplyr select inner_join rename_with
#' @importFrom reshape2 melt dcast

get_vaf_df = function(vaf_file, obj, min_ccf=0) {
  lineages = obj$lineages
  highlight = select_relevant_clusters(obj, min_ccf=min_ccf)

  dataframe = obj$dataframe %>% filter(labels %in% highlight) %>%
    reshape2::melt() %>% tidyr::separate(variable, into=c("else", "timepoint", "lineage"), sep="_") %>%
    mutate("else"=NULL) %>%
    tidyr::pivot_wider(names_from="timepoint", values_from="value", names_prefix="cov_")
  IS_keep = dataframe$IS

  vaf_ref_alt = read.csv(vaf_file) %>% filter(IS %in% IS_keep, lineage %in% lineages)

  vaf_cov_df = dplyr::inner_join(vaf_ref_alt, dataframe, by=c("IS", "lineage")) %>%
    tidyr::pivot_longer(cols=starts_with("dp.ref.alt"), names_to="timepoint", values_to="dp.ref.alt") %>%
    tidyr::separate(timepoint, into=c("else", "timepoint"), sep="_") %>%
    tidyr::pivot_longer(cols=starts_with("vaf"), names_to="timepoint2", values_to="vaf") %>%
    tidyr::separate(timepoint2, into=c("else", "timepoint2"), sep="_") %>%
    filter(timepoint==timepoint2) %>% mutate("else"=NULL, timepoint2=NULL) %>%
    dplyr::select(-starts_with("dp_"), -starts_with("cov_")) %>%
    tidyr::separate("dp.ref.alt", into=c("dp", "ref", "alt"), sep="[:]") %>%
    mutate(dp=as.integer(dp), ref=as.integer(ref), alt=as.integer(alt)) %>%
    mutate(dp=ifelse(dp<alt, alt+ref, dp))

  return(vaf_cov_df)
}


# Function to get from a vaf dataframe obtained by get_vaf_df() the input for a VIBER run
get_input_viber = function(vaf_df) {
  trials = vaf_df %>% reshape2::dcast(IS+labels+mutation+experiment ~ timepoint+lineage, value.var="dp")
  successes = vaf_df %>% reshape2::dcast(IS+labels+mutation+experiment ~ timepoint+lineage, value.var="alt")

  joined = dplyr::inner_join(trials, successes, by=c("mutation","IS","labels","experiment"),
                      suffix=c(".dp", ".alt"))
  return(list("successes"=successes, "trials"=trials, "joined"=joined))
}


# As input a mvnmm object with already the VAF dataframe
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


