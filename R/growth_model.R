fit_growth = function(x,
                      steps=500,
                      highlight=list(),
                      timepoints_to_int=list(),
                      growth="",
                      force=T) {

  if (purrr::is_empty(highlight)) highlight = x %>% get_highlight(highlight=highlight, mutations=F)
  timepoints_to_int = map_timepoints_int(x, timepoints_to_int)

  pop_df = x %>%
    get_muller_pop(mutations=T, timepoints_to_int=timepoints_to_int)

  rates.df = x %>% get_growth_rates() # data.frame()
  evaluated = c()
  if (!purrr::is_empty(rates.df)) evaluated = rates.df$Identity %>% unique()

  if (force & !purrr::is_empty(rates.df)) {
    rates.df = rates.df %>% filter(!Identity %in% highlight)
    evaluated = rates.df$Identity %>% unique()
  }

  for (cluster in highlight) {

    if (!cluster %in% evaluated) {

      pop_df.cl = pop_df %>%
        filter(grepl(paste(cluster,".",sep=""), Identity) | Identity == cluster) %>%
        dplyr::select(-starts_with("lm"))
      parents = get_parents(x, highlight=cluster)

      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=cluster,
                                   pop_df=pop_df.cl,
                                   parents=parents,
                                   growth=growth,
                                   steps=steps,
                                   clonal=TRUE,
                                   timepoints_to_int=timepoints_to_int)

      # fit the growth rate in all the parents
      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=parents$Parent %>% unique(),
                                   pop_df=pop_df,
                                   parents=parents,
                                   growth=growth,
                                   steps=steps,
                                   clonal=F,
                                   timepoints_to_int=timepoints_to_int)
      # then in all the children
      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=parents$Identity %>% unique(),
                                   pop_df=pop_df,
                                   parents=parents,
                                   growth=growth,
                                   steps=steps,
                                   clonal=F,
                                   timepoints_to_int=timepoints_to_int)

      evaluated = c(cluster, evaluated)
    }
  }

  x$growth_rates = rates.df
  x$tp_to_int = timepoints_to_int

  return(x)
}


get_growth_rates = function(x) {
  if ("growth_rates" %in% names(x))
    return(x$growth_rates)
  return(data.frame())
}


fit_growth_multiple_clones = function(rates.df,
                                      clusters,
                                      pop_df,
                                      parents,
                                      growth,
                                      steps,
                                      clonal,
                                      timepoints_to_int) {
  for (subcl in clusters) {

    if (!subcl %in% (rates.df$Identity)) {

      if (purrr::is_empty(rates.df))
        rates.df = fit_growth_single_clone(pop_df=pop_df,
                                           cluster=subcl,
                                           timepoints_to_int=timepoints_to_int,
                                           p.rates=get_parent_rate(parents, rates.df, subcl),
                                           clonal=clonal,
                                           growth=growth,
                                           steps=steps)
      else rates.df = rates.df %>%
          dplyr::add_row(
            fit_growth_single_clone(pop_df=pop_df,
                                    cluster=subcl,
                                    timepoints_to_int=timepoints_to_int,
                                    p.rates=get_parent_rate(parents, rates.df, subcl),
                                    clonal=clonal,
                                    growth=growth,
                                    steps=steps)
          )
    }
  }
  return(rates.df)
}


fit_growth_single_clone = function(pop_df,
                                   cluster,
                                   timepoints_to_int,
                                   steps=500,
                                   p.rates=list("exp"=NULL, "log"=NULL),
                                   clonal=FALSE,
                                   growth="",
                                   random_state=25) {
  torch = reticulate::import("torch")
  py_pkg = reticulate::import("pylineaGT")
  # reticulate::source_python("../pyLineaGT/pylineaGT/explogreg.py")

  input.clone = pop_df %>%
    filter(Identity == cluster) %>%
    dplyr::select(-Frequency) %>%
    dplyr::mutate(Population=ifelse((Generation == 0 & clonal), 1, Population)) %>%
    pivot_wider(id_cols="Generation", names_from="Lineage", values_from="Population")

  times = input.clone$Generation; times = torch$tensor(times)$int()$unsqueeze(as.integer(1))
  y = input.clone[2:ncol(input.clone)] %>% as.matrix(); y = torch$tensor(y)

  x.reg = py_pkg$explogreg$Regression(times, y)

  p.rate.exp = p.rate.log = NULL
  ## exp training
  if (growth=="" | growth=="exp") {
    if (!is.null(p.rates[["rate.exp"]])) p.rate.exp = torch$tensor(p.rates[["rate.exp"]])$float()

    losses.exp = x.reg$train(regr="exp", p_rate=p.rate.exp, steps=as.integer(steps), random_state=as.integer(random_state))
    p.exp = x.reg$get_learned_params()
  }

  ## log training
  if (growth=="" | growth=="log") {
    if (!is.null(p.rates[["rate.log"]])) p.rate.log = torch$tensor(p.rates[["rate.log"]])$float()

    losses.log = x.reg$train(regr="log", p_rate=p.rate.log, steps=as.integer(steps), random_state=as.integer(random_state))
    p.log = x.reg$get_learned_params()
  }

  return( get_growth_params(timepoints_to_int, p.exp, p.log, lineages=pop_df$Lineage %>% unique(), cluster=cluster) )
}


get_parent_rate = function(parents, rates.df, cluster) {
  if (!cluster %in% parents$Identity)
    return(list("exp"=NULL, "log"=NULL))

  par = parents %>% filter(Identity==cluster) %>% dplyr::pull(Parent)  ## parent of subcl
  par.rates = rates.df %>%
    filter(Identity==par) %>%
    dplyr::select(dplyr::contains("rate.exp"), dplyr::contains("rate.log"), -dplyr::contains("p_rate")) %>%
    as.list()

  if (NA %in% unlist(par.rates))
    return(list("exp"=NULL, "log"=NULL))
  return(par.rates)
}


get_growth_params = function(timepoints_to_int,
                             rates.exp=NULL,
                             rates.log=NULL,
                             lineages=list(),
                             cluster=NULL) {

  params.exp = params.log = data.frame()

  if (!is.null(rates.exp))
    params.exp = params.exp %>%
      tibble::add_column("fitness.exp"=NA, "init_t.exp"=NA, "Lineage"=as.character(NA), "p_rate.exp"=NA, "sigma.exp"=NA) %>%
      tibble::add_row("Lineage"=lineages,
                      "p_rate.exp"=rates.exp$parent_rate,
                      "fitness.exp"=rates.exp$fitness,
                      "init_t.exp"=rates.exp$init_time) %>%
      dplyr::mutate("rate.exp"=ifelse(is.na(p_rate.exp),
                                      fitness.exp,
                                      p_rate.exp * (1+fitness.exp)),
                    "Identity"=cluster,
                    "sigma.exp"=list(setNames(object=rates.exp$sigma, nm=c(0, unlist(timepoints_to_int))))) %>%
      tibble::as_tibble()

  print(rates.log)
  if (!is.null(rates.log))
    params.log = params.log %>%
      tibble::add_column("fitness.log"=NA, "K.log"=NA, "init_t.log"=NA, "Lineage"=as.character(NA),
                         "p_rate.log"=NA, "sigma.log"=NA) %>%
      tibble::add_row("fitness.log"=rates.log$fitness,
                      "K.log"=rates.log$carr_capac,
                      "init_t.log"=rates.log$init_time,
                      "Lineage"=lineages,
                      "p_rate.log"=rates.log$parent_rate) %>%
      dplyr::mutate("rate.log"=ifelse(is.na(p_rate.log),
                                      fitness.log,
                                      p_rate.log * (1+fitness.log)),
                    "Identity"=cluster,
                    "sigma.log"=list(setNames(object=rates.log$sigma, nm=c(0, unlist(timepoints_to_int))))) %>%
      tibble::as_tibble()

  if (is.null(rates.exp))
    return(params.log)

  if (is.null(rates.log))
    return(params.exp)

  return(
    dplyr::inner_join(params.exp, params.log, by=c("Lineage", "Identity"))
  )
}


