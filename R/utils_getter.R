get_mean_long = function(x) {
  return(
    x %>%
      get_mean() %>%
      as.data.frame() %>% tibble::rownames_to_column(var="labels") %>%
      dplyr::mutate(labels=factor(labels, levels=unique(labels))) %>%
      tidyr::pivot_longer(cols=starts_with("cov"), names_to="cov.timepoints.lineage", values_to="mean_cov") %>%
      tidyr::separate(cov.timepoints.lineage, into=c("else", "timepoints", "lineage"), sep="[.]") %>%
      dplyr::mutate("else"=NULL) %>%
      mutate_tp(fn=as.integer, colnm="timepoints")
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
      dplyr::mutate("else"=NULL) %>%
      mutate_tp(fn=as.integer, colnm="timepoints")
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
