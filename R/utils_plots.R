get_pairs = function(dataset, columns) {
  comb = t(combn(names(dataset %>% dplyr::select(dplyr::all_of(columns))), 2)) %>%
    as.data.frame(stringsAsFactors=FALSE) %>%
    mutate(pair_name=paste(V1, V2, sep=":"))
  return(comb)
}


# To compute the Gaussian multivariate density given a fitted object
compute_density = function(x) {
  mean = get_mean(x)
  sigma = get_covariance_Sigma(x)

  density = data.frame()
  for (cl in get_unique_labels(x)) {
    dd = as.data.frame(MASS::mvrnorm(n=1000, mu=mean[cl,], Sigma=sigma[[cl]]))
    colnames(dd) = x$dimensions
    dd$labels = cl
    density = rbind(density, dd)
  }
  density = density %>% mutate(labels=factor(labels, levels=get_unique_labels(x)))
  return(density)
}


split_to_camelcase = function(txt) {
  txt = stringr::str_replace_all(txt, "\\_|\\.", " ")
  return(paste(toupper(substring(txt,1,1)), substring(txt,2), sep=""))
}


get_colors = function(x=NULL, list_lab=list(), color_palette=list()) {
  if (purrr::is_empty(list_lab)) {
    N = x$K
    colss = Polychrome::createPalette(N, c("#856de3", "#9e461c"), target="normal", range=c(15, 80), M=1000)
    colss = colss[1:N]
    try({ names(colss) = x$params$labels %>% levels() }, silent=T) }
  else {
    # means we want colors for the subclones
    colss = c()
    for (cl in color_palette %>% names) {
      mut_cl = get_unique_muts_labels(x, clusters=c(cl))
      if (!purrr::is_empty(mut_cl)) {
        n_cols = mut_cl %>% length() + 1
        new_cols = Polychrome::createPalette(n_cols, c(color_palette[cl]),
                                             target="normal", range=c(15, 80), M=1000)[2:n_cols]
        names(new_cols) = mut_cl
        colss = c(colss, new_cols)
      }
    }
    # N = list_lab %>% length()
    # colss = Polychrome::createPalette(N, c("#856de3", "#9e461c"), target="normal", range=c(15, 80), M=100000)
    # colss = colss[1:N]
    # names(colss) = list_lab
  }
  return(colss)
}


highlight_palette = function(color_palette, highlight=c()) {
  if (purrr::is_empty(highlight)) return(color_palette)

  remove = color_palette[!names(color_palette)%in% highlight] %>% names
  keep = color_palette[names(color_palette)%in% highlight]
  grey_col = gray(runif(remove %>% length(), 0.6, 0.8))
  names(grey_col) = remove
  return(c(keep, grey_col))
}


select_relevant_clusters = function(x, min_frac) {
  return(
    x %>%
      get_muller_pop() %>%
      group_by(Identity) %>%
      filter(any(Frequency > min_frac), Identity!="P") %>%
      dplyr::pull(Identity) %>%
      unique()
  )
}


retrieve_clusters = function(x, min_frac, highlight) {
  if (purrr::is_empty(highlight)) highlight = x %>% get_unique_labels()
  return(intersect(select_relevant_clusters(x, min_frac), highlight))
}



# reshape_vaf_dataframe_long = function(x) {
#   vaf = x %>% get_vaf_dataframe() %>% mutate(labels_mut=paste(labels,labels_viber,sep=".")) %>%
#     dplyr::select(starts_with("vaf"), mutation, IS, contains("labels"), contains("viber")) %>%
#     tidyr::pivot_longer(cols=starts_with("vaf"), names_to="timepoints_lineage", values_to="vaf") %>%
#     separate(timepoints_lineage, into=c("vv","timepoints","lineage")) %>%
#     mutate(timepoints=paste(vv,timepoints,sep="."),vv=NULL) %>%
#     tidyr::pivot_wider(names_from=timepoints, values_from="vaf")
#
#   try(expr = {vaf = vaf %>% dplyr::select(-"vaf.over")}, silent=T)
#   try(expr = {vaf = vaf %>% dplyr::select(-"vaf.steady")}, silent=T)
#
#   return(vaf)
# }

