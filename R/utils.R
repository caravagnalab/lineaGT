update_color_palette = function(x, clusters=c(), label="") {
  list_lab = x %>%
    get_unique_muts_labels(clusters=clusters, label=label)
  return(
    c(get_color_palette(x, label)[x %>% get_unique_labels()],
      get_colors(x=x,
                 list_lab=list_lab,
                 color_palette=get_color_palette(x, label)[x %>% get_unique_labels()],
                 label=label))
  )
}


get_highlight = function(x, min_frac=0, highlight=c(), mutations=F, label="") {
  if (mutations) {
    if (purrr::is_empty(highlight))
      highlight = select_relevant_clusters(x, min_frac)
    highlight_v = get_unique_muts_labels(x, highlight, label=label)

    return(
      c(c(highlight, highlight_v))
    )
  }

  if (purrr::is_empty(highlight)) highlight = x %>% get_unique_labels()
  return(
    intersect(select_relevant_clusters(x, min_frac), highlight)
  )
}


select_relevant_clusters = function(x, min_frac) {
  return(
    x %>%
      get_muller_pop(timepoints_to_int=NULL, exp_coef=F) %>%
      group_by(Identity) %>%
      filter(any(Frequency > min_frac), Identity!="P") %>%
      dplyr::pull(Identity) %>%
      unique()
  )
}

