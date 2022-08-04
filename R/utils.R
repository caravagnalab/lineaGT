update_color_palette = function(x, clusters=c()) {
  list_lab = x %>%
    get_unique_muts_labels(clusters=clusters)

  color_palette = x %>% get_color_palette()
  unq = get_unique_labels(x)

  return(
    c(color_palette[unq],
      get_colors(x=x,
                 list_lab=list_lab,
                 color_palette=color_palette[unq]))
  )
}


get_highlight = function(x, min_frac=0, highlight=c(), mutations=F) {
  if (mutations) {
    if (purrr::is_empty(highlight)) highlight = select_relevant_clusters(x, min_frac)
    highlight_v = get_unique_muts_labels(x, highlight)

    return( c(highlight, highlight_v) )
  }

  if (purrr::is_empty(highlight)) highlight = x %>% get_unique_labels()
  return( intersect(select_relevant_clusters(x, min_frac), highlight) )
}


select_relevant_clusters = function(x, min_frac) {
  return(
    x %>%
      get_muller_pop(force=F) %>%
      dplyr::group_by(Identity) %>%
      dplyr::filter(any(Frequency > min_frac), Identity!="P") %>%
      dplyr::pull(Identity) %>%
      unique()
  )
}

