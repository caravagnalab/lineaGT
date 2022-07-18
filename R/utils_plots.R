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


get_colors = function(x=NULL, list_lab=list(), color_palette=list(), label="") {
  if (purrr::is_empty(list_lab)) {
    N = x$K
    colss = Polychrome::createPalette(N, c("#856de3", "#9e461c"), target="normal", range=c(15, 80), M=1000)
    colss = colss[1:N]
    try({ names(colss) = x$params$labels %>% levels() }, silent=T) }
  else {
    # means we want colors for the subclones
    colss = c()
    for (cl in color_palette %>% names) {
      mut_cl = get_unique_muts_labels(x, clusters=c(cl), label=label)
      if (!purrr::is_empty(mut_cl)) {
        n_cols = mut_cl %>% length() + 1
        new_cols = Polychrome::createPalette(n_cols, c(color_palette[cl]),
                                             target="normal", range=c(15, 80), M=1000)[2:n_cols]
        names(new_cols) = mut_cl
        colss = c(colss, new_cols)
      }
    }
  }
  return(colss)
}


highlight_palette = function(x, highlight=c(), label="") {
  if (purrr::is_empty(highlight)) return(x$color_palette)
  color_palette = get_color_palette(x, label)

  remove = color_palette[!names(color_palette)%in% highlight] %>% names
  keep = color_palette[names(color_palette)%in% highlight]
  grey_col = gray(runif(remove %>% length(), 0.6, 0.8))
  names(grey_col) = remove
  return(c(keep, grey_col))
}


map_timepoints_int = function(x, timepoints_to_int=list()) {
  if (!purrr::is_empty(timepoints_to_int)) return(timepoints_to_int)

  if (!purrr::is_empty(x %>% get_tp_to_int())) return(x %>% get_tp_to_int())

  tp = x %>% get_timepoints()
  names(tp) = as.character(tp)

  # if is numeric or integer
  if (is.numeric(tp)) {
    tp = tp %>% sort() %>% as.list()
    x$tp_to_int = tp
    return(tp)
  }
  else {
    tryCatch(expr =
               {
                 tp = as.numeric(tp) %>% sort() %>% as.list()
                 names(tp) = as.character(tp)
                 x$tp_to_int = tp
                 return(tp)
               },
             warning = function(w) {} )
  }

  if (is.factor(tp)) {
    names(tp) = tp %>% levels()
    cli::cli_warn("The provided timepoints are Factors.
                  They will be converted to integer, with time unit of 50.")
    tp = seq(from=50, to=50*length(tp), length.out=length(tp)) %>% as.list()
    x$tp_to_int = tp
    return(tp)
  }

  if (is.character(tp)) {
    cli::cli_warn("The provided timepoints are character.
                  If you want to provide a temporal order, insert them as numeric or factors.
                  They will be converted to integer, with time unit of 50.")
    tp = seq(from=50, to=50*length(tp), length.out=length(tp)) %>% as.list()
    names(tp) = x %>% get_timepoints()
    x$tp_to_int = tp
    return(tp)
  }
}



