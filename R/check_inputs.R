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

