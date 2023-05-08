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
  if (mutations && have_vaf_df(x)) {
    if (purrr::is_empty(highlight)) highlight = select_relevant_clusters(x, min_frac)
    highlight_v = get_unique_muts_labels(x, highlight)

    return( unique(c(highlight, highlight_v)) )
  }

  if (purrr::is_empty(highlight)) highlight = x %>% get_unique_labels()
  return( intersect(select_relevant_clusters(x, min_frac), highlight) )
}


select_relevant_clusters = function(x, min_frac) {
  return(
    x %>%
      get_mean_long() %>%
      dplyr::group_by(timepoints, lineage) %>%
      dplyr::mutate(frac=mean_cov/sum(mean_cov)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(labels) %>%
      dplyr::filter(any(frac > min_frac), labels!="P") %>%
      dplyr::pull(labels) %>%
      unique() %>%
      as.character()
  )
}


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
    dd = as.data.frame(MASS::mvrnorm(n=500, mu=mean[cl,], Sigma=sigma[[cl]], empirical=T))
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
    colss = Polychrome::createPalette(N, c("#856de3","#9e461c"), target="normal", range=c(15,80), M=1000)[1:N]
    try(expr = { names(colss) = get_unique_labels(x) }, silent=T)
    return(colss)
  }

  # means we want colors for the subclones
  colss = c()
  for (cl in names(color_palette)) {
    mut_cl = get_unique_muts_labels(x, clusters=cl)
    if (!purrr::is_empty(mut_cl)) {
      n_cols = length(mut_cl) + 1
      new_cols = Polychrome::createPalette(n_cols, c(color_palette[cl]),
                                           target="normal", range=c(15, 80), M=1000)[2:n_cols] %>%
        setNames(nm=mut_cl)
      colss = c(colss, new_cols)
    }
  }
  return(colss)
}


highlight_palette = function(x, highlight=c()) {
  if (purrr::is_empty(highlight)) return(x$color_palette)
  color_palette = get_color_palette(x)

  remove = color_palette[!names(color_palette)%in% highlight] %>% names
  keep = color_palette[names(color_palette)%in% highlight]
  grey_col = gray(runif(remove %>% length(), 0.6, 0.8))
  names(grey_col) = remove
  return(c(keep, grey_col))
}


map_timepoints_int = function(x, timepoints_to_int=list()) {
  if (!purrr::is_empty(timepoints_to_int)) return(timepoints_to_int %>% unlist())

  suppressMessages( expr = { if (!purrr::is_empty(x %>% get_tp_to_int())) return(x %>% get_tp_to_int()) } )

  tp = x %>% get_timepoints()
  names(tp) = as.character(tp)

  # if is numeric or integer
  if (is.numeric(tp)) return( tp %>% sort() )
  else # check if they are convertible to numeric
    tryCatch(expr = {
      tp = as.numeric(tp) %>% sort()
      return( tp %>% setNames( as.character(tp) ) )
      },
      warning = function(w) {} )

  if ( is.factor(tp) ) {
    cli::cli_alert_warning("The provided timepoints are Factors.
                            They will be converted to integer, with time unit of 50.")
    tp.int = seq(from=50, to=50*length(tp), length.out=length(tp)) %>%
      setNames( nm=levels(tp) )
    return( tp.int )
  }

  if (is.character( tp )) {
    cli::cli_alert_warning("The provided timepoints are characters.
                            If you want to provide a temporal order, insert them as numeric or factors.
                            They will be converted to integer, with time unit of 50.")
    tp.int = seq(from=50, to=50*length(tp), length.out=length(tp)) %>%
      setNames( nm=tp )
    return( tp.int )
  }
}


mutate_tp = function(dataset, fn, colnm="timepoints") {
  tryCatch(expr = {
    dataset = dataset %>%
      dplyr::mutate(dplyr::across(!!colnm, fn))
    return(dataset)
  },
  warning = function(w) return(dataset),
  error = function(e) return(dataset)
  )

  return(dataset)
}



