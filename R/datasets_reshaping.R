long_to_wide_cov = function(cov.df) {
  ## transforms the input dataset from long to wide format
  ## input columns are: "coverage", "timepoints", "lineage", "IS"

  return(
    cov.df %>%
      tidyr::pivot_wider(names_from=c("timepoints","lineage"),
                         names_prefix="cov.",
                         names_sep=".",
                         values_from="coverage",
                         values_fill=0)  # fill missing values with 0 -> means 0 coverage
    )
}


long_to_wide_muts = function(vaf.df) {
  return(
    vaf.df %>%
      dplyr::select(timepoints, lineage, IS, mutation,
                    dplyr::starts_with("alt"),
                    dplyr::starts_with("ref"),
                    dplyr::starts_with("dp"),
                    dplyr::starts_with("vaf"),
                    dplyr::starts_with("theta"),
                    dplyr::starts_with("labels"),
                    dplyr::starts_with("pi")) %>%

      # keep "IS" and "mutation" as id columns
      tidyr::pivot_wider(names_from=c("timepoints","lineage"), names_sep=".",
                         values_from=c(dplyr::starts_with("alt"),
                                       dplyr::starts_with("ref"),
                                       dplyr::starts_with("dp"),
                                       dplyr::starts_with("vaf"),
                                       dplyr::starts_with("theta")),
                         values_fill=0)  # fill missing values with 0 -> means 0 ref/alt/vaf
  )
}


wide_to_long_cov = function(cov.df) {
  return(
    cov.df %>%
      tidyr::pivot_longer(cols=starts_with("cov"),
                          names_to="else.time.lineage",
                          values_to="coverage") %>%
      tidyr::separate("else.time.lineage",
                      into=c("else","timepoints","lineage"),
                      sep="[.]") %>%
      dplyr::mutate("else"=NULL)
  )
}


wide_to_long_muts = function(vaf.df) {
  return(
    vaf.df %>%
      tidyr::pivot_longer(cols=c(dplyr::starts_with("alt"),
                                 dplyr::starts_with("ref"),
                                 dplyr::starts_with("dp"),
                                 dplyr::starts_with("theta"),
                                 dplyr::starts_with("vaf")),
                          names_to="type.timepoints.lineage") %>%
      tidyr::separate("type.timepoints.lineage",
                      into=c("type", "timepoints", "lineage"),
                      sep="[.]") %>%
      tidyr::pivot_wider(names_from="type",
                         values_from="value")
  )
}


