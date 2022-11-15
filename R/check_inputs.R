have_muts_fit = function(x) {
  x.muts = suppressMessages(x %>% get_muts_fit())
  if (purrr::is_empty(x.muts))
    return(FALSE)
  return(TRUE)
}


have_vaf_df = function(x) {
  vaf.df = suppressMessages(x %>% get_vaf_dataframe())
  if (purrr::is_empty(vaf.df))
    return(FALSE)
  return(TRUE)
}


have_growth_rates = function(x) {
  rates.df = suppressMessages(x %>% get_growth_rates())
  if (purrr::is_empty(rates.df))
    return(FALSE)
  return(TRUE)
}


have_phylo_fits = function(x) {
  trees = suppressMessages(x %>% get_trees())
  if (purrr::is_empty(trees))
    return(FALSE)
  return(TRUE)
}


have_tp_int = function(x) {
  tp.int = suppressMessages(x %>% get_tp_to_int())
  if (purrr::is_empty(tp.int))
    return(FALSE)
  return(TRUE)
}


have_pop_df = function(x) {
  if ("population.df" %in% names(x))
    return(TRUE)
  return(FALSE)
}


have_pop_df_muts = function(x) {
  muts_labs = get_unique_muts_labels(x)
  if ("population.df" %in% names(x) && sum(muts_labs %in% get_pop_df(x)$Identity)==length(muts_labs))
    return(TRUE)
  return(FALSE)
}


have_corrected_pops = function(x, estimate_npops) {
  if (have_pop_df(x) && !estimate_npops)
    return(TRUE)

  if (have_pop_df(x) && estimate_npops && ("Population.corr" %in% colnames(get_pop_df(x))))
    return(TRUE)

  return(FALSE)

}

