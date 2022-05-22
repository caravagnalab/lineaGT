
long_to_wide_cov = function(dataset) {
  ## transforms the input dataset from long to wide format
  ## input columns are: "coverage", "timepoints", "lineage", "IS"

  return(
    dataset %>%
      tidyr::pivot_wider(names_from=c("timepoints","lineage"), names_prefix="cov.", names_sep=".",
                         values_from="coverage", values_fill=0)
    )
}


long_to_wide_muts = function(vaf.df) {
  return(
    vaf.df %>%
      dplyr::select(alt, dp, vaf, timepoints, lineage, IS, mutation,
                    dplyr::starts_with("labels"), dplyr::contains("pi")) %>%
      tidyr::pivot_wider(names_from=c("timepoints","lineage"), names_sep=".",
                         values_from=c("alt","dp","vaf"))
  )
}


wide_to_long_cov = function(dataset) {
  return(
    dataset %>%
      tidyr::pivot_longer(cols=starts_with("cov"), names_to="else.time.lineage", values_to="coverage") %>%
      separate(else.time.lineage, into=c("else","timepoints","lineage")) %>%
      mutate("else"=NULL)
  )
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

