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


get_colors = function(x=NULL, list_lab=list()) {
  if (purrr::is_empty(list_lab)) {
    N = x$K
    colss = Polychrome::createPalette(N, c("#856de3", "#9e461c"), target="normal", range=c(15, 80), M=100000)
    colss = colss[1:N]
    try({ names(colss) = x$params$labels %>% levels() }, silent=T) }
  else {
    N = list_lab %>% length()
    colss = Polychrome::createPalette(N, c("#856de3", "#9e461c"), target="normal", range=c(15, 80), M=100000)
    colss = colss[1:N]
    names(colss) = list_lab
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



get_muller_pop = function(x, means=list()) {
  if (purrr::is_empty(means)) means = get_mean(x)

  pop_df = means %>% as.data.frame() %>%
    tibble::rownames_to_column(var="Identity") %>%
    tidyr::pivot_longer(cols=c(starts_with("cov"), starts_with("vaf")), names_to="timepoints_lineage", values_to="Population") %>%
    tidyr::separate(timepoints_lineage, into=c("else", "Generation", "Lineage"), sep="\\.|\\_") %>%
    mutate("else"=NULL, Population=ifelse(Population==0, 0.001, Population)) %>%
    group_by(Generation, Lineage) %>%
    mutate(Frequency=Population/sum(Population)) %>%
    dplyr::ungroup()

  pop_df = rbind(pop_df, list("Identity"=rep("P", x$`T`+x$lineages %>% length()),
                              "Generation"=c(x$dimensions,rep("init", x$lineages %>% length())),
                              "Population"=rep(1, x$`T`+x$lineages %>% length()),
                              "Lineage"=rep(x$lineages, 3+1),
                              "Frequency"=rep(1, x$`T`+x$lineages %>% length())))
  pop_df = rbind(pop_df, list("Identity"=rep(x %>% get_unique_labels(), x$lineages %>% length()),
                       "Generation"=rep("init", x$K*(x$lineages %>% length())),
                       "Population"=rep(1, x$K*(x$lineages %>% length())),
                       "Lineage"=rep(x$lineages, each=x$K),
                       "Frequency"=rep(1/x$K, x$K*(x$lineages %>% length())))) %>%
    mutate(Generation=dplyr::case_when(grepl("early", Generation) ~ "60",
                                       grepl("mid", Generation) ~ "140",
                                       grepl("late", Generation) ~ "280",
                                       grepl("init", Generation) ~ "1")) %>%
    mutate(Generation=as.numeric(Generation)) %>%
    group_by(Identity, Lineage) %>%
    mutate(lm_a=coef(lm(log1p(Population)~Generation))[1],
           lm_r=coef(lm(log1p(Population)~Generation))[2]) %>% ungroup()
  return(pop_df)
}


get_muller_edges = function(x, labels=list()) {
  if (purrr::is_empty(labels)) return(data.frame("Parent"="P", "Identity"=get_unique_labels(x)))
  return(data.frame("Parent"="P", "Identity"=labels))
}


select_relevant_clusters = function(x, min_frac, means=list()) {
  pop_df = get_muller_pop(x, means)
  clusters_keep = (pop_df %>% group_by(Identity) %>%
                     filter(any(Frequency > min_frac), Identity!="P"))$Identity %>% unique()
  return(clusters_keep)
}


reshape_vaf_dataframe_long = function(x) {
  vaf = x %>% get_vaf_dataframe() %>% mutate(labels_mut=paste(labels,labels_viber,sep=".")) %>%
    dplyr::select(starts_with("vaf"), mutation, IS, contains("labels"), contains("viber")) %>%
    tidyr::pivot_longer(cols=starts_with("vaf"), names_to="timepoints_lineage", values_to="vaf") %>%
    separate(timepoints_lineage, into=c("vv","timepoints","lineage")) %>%
    mutate(timepoints=paste(vv,timepoints,sep="."),vv=NULL) %>%
    tidyr::pivot_wider(names_from=timepoints, values_from="vaf")

  try(expr = {vaf = vaf %>% dplyr::select(-"vaf.over")}, silent=T)
  try(expr = {vaf = vaf %>% dplyr::select(-"vaf.steady")}, silent=T)

  return(vaf)
}

# reshape_theta_long = function(x) {
#   theta = get_binomial_theta(x) %>%
#     tidyr::pivot_longer(cols=starts_with("vaf"), names_to="timepoints_lineage", values_to="vaf") %>%
#     tidyr::separate(timepoints_lineage, into=c("timepoints","lineage"), sep="_") %>%
#     tidyr::pivot_wider(names_from="timepoints", values_from="vaf") %>%
#     mutate(labels=labels_mut) %>% separate(labels, into=c("labels","else"), sep="[.]") %>%
#     mutate("else"=NULL)
#
#   return(theta)
# }
#
