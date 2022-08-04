#' Infer growth rates for each clone and subclone.
#'
#' @param x a mvnmm object.
#' @param steps maximum number of steps for inference.
#' @param highlight set of clusters to run the inference for.
#' If not specified, it will be run on all the clusters.
#' @param timepoints_to_int if the provided timepoints are not integers nor a timepoints-to-int list is
#' stored in \code{x}, a list mapping their values to integers is required.
#' @param growth_model string specifying the type of growth model to use, between \code{exp} and \code{log}
#' corresponding to exponential and logistic models, respectively.
#' @param force if the model has already been fitted, setting \code{force} to \code{FALSE} will keep the
#' computed rates. Setting \code{force} to \code{TRUE} will fit the model again for the specified clusters.
#'
#' @return
#'
#' @importFrom purrr is_empty
#' @importFrom dplyr select filter add_row mutate pull contains select
#'
#' @export

fit_growth_rates = function(x,
                            steps=500,
                            highlight=c(),
                            timepoints_to_int=c(),
                            growth_model="",
                            force=T,
                            tree_score=1) {

  highlight = get_highlight(x, highlight=highlight, mutations=F)
  timepoints_to_int = map_timepoints_int(x, timepoints_to_int)

  if (have_muts_fit(x)) mutations = T else mutations = F

  pop_df = x %>%
    get_muller_pop(mutations=mutations, timepoints_to_int=timepoints_to_int, tree_score=tree_score)

  rates.df = x %>% get_growth_rates()
  if (have_growth_rates(x)) evaluated = rates.df$Identity %>% unique() else evaluated = list()

  # if force is TRUE, then all the clusters are fitted again
  if (force & have_growth_rates(x)) {
    rates.df = rates.df %>% filter(!Identity %in% highlight)
    evaluated = rates.df$Identity %>% unique()
  }

  for (cluster in highlight) {
    if (!cluster %in% evaluated) {

      cli::cli_process_start(paste0("Starting growth models inference of clone ", cluster))

      # filter the dataset
      pop_df.cl = pop_df %>%
        dplyr::filter(grepl(paste(cluster,".",sep=""), Identity) | Identity == cluster) %>%
        dplyr::select(-starts_with("lm"))
      # get the identity-parents dataframe of clone "cluster"
      parents = get_parents(x, highlight=cluster)

      # first in the clonal cluster of ISs
      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=cluster,
                                   pop_df=pop_df.cl,
                                   parents=parents,
                                   growth=growth_model,
                                   steps=steps,
                                   clonal=TRUE,
                                   timepoints_to_int=timepoints_to_int)

      # fit the growth rate in all the parents
      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=parents$Parent %>% unique(),
                                   pop_df=pop_df,
                                   parents=parents,
                                   growth=growth_model,
                                   steps=steps,
                                   clonal=F,
                                   timepoints_to_int=timepoints_to_int)

      # then in all the children
      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=parents$Identity %>% unique(),
                                   pop_df=pop_df,
                                   parents=parents,
                                   growth=growth_model,
                                   steps=steps,
                                   clonal=F,
                                   timepoints_to_int=timepoints_to_int)

      cli::cli_process_done()

      evaluated = c(cluster, evaluated)
    }
  }

  x = add_growth_rates(x, rates.df)
  x = add_tp_int(x, timepoints_to_int)

  return(x)
}


sort_clusters_edges = function(clusters, parents) {
  # clusters is a list of clusters
  # parents is the edges tibble
  nn = c(parents$Parent, parents$Identity) %>% unique() %>% length()
  nodes = c(parents$Parent, parents$Identity) %>% unique()

  root = setdiff(parents$Parent, parents$Identity)
  node.name = parents %>% dplyr::filter(Parent==root) %>% dplyr::pull(Identity)
  cls.sort = c(root)

  for (i in 1:nn) {
    if (node.name %in% clusters) cls.sort = c(cls.sort, node.name)
    node.name = parents %>% dplyr::filter(Identity==node.name) %>% dplyr::pull(Identity)
  }

  return(
    intersect(cls.sort, clusters) %>% unique()
  )
}


fit_growth_multiple_clones = function(rates.df,
                                      clusters,
                                      pop_df,
                                      parents,
                                      growth,
                                      steps,
                                      clonal,
                                      timepoints_to_int) {

  if (!clonal & nrow(parents)==0) return(rates.df)

  # get the clusters sorted by tree hierarchy
  if (!clonal) clusters = sort_clusters_edges(clusters, parents)

  for (subcl in clusters) {

    if (!subcl %in% (rates.df$Identity) && (subcl %in% pop_df$Identity)) {
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
                                   random_state=25,
                                   py_pkg=NULL) {

  torch = reticulate::import("torch")
  if (is.null(py_pkg)) py_pkg = reticulate::import("pylineaGT")

  input.clone = pop_df %>%
    filter(Identity == cluster) %>%
    dplyr::select(-Frequency) %>%
    dplyr::mutate(Population=ifelse((Generation == 0 & clonal), 1, Population)) %>%
    tidyr::pivot_wider(id_cols="Generation", names_from="Lineage", values_from="Population")

  times = input.clone$Generation
  y = input.clone[2:ncol(input.clone)]
  lineages = colnames(y)

  times = torch$tensor(times)$int()$unsqueeze(as.integer(1))
  y = torch$tensor(y %>% as.matrix())

  p.rate.exp = p.rate.log = NULL

  x.reg = py_pkg$explogreg$Regression(times, y)
  if (growth=="" | growth=="exp") {  # exp training
    if (!is.null(p.rates[["rate.exp"]]))
      p.rate.exp = torch$tensor(p.rates[["rate.exp"]])$float()

    losses.exp = x.reg$train(regr="exp", p_rate=p.rate.exp, steps=as.integer(steps), random_state=as.integer(random_state))
    p.exp = x.reg$get_learned_params()
    ll.exp = x.reg$compute_log_likelihood() %>% setNames(nm=lineages)
  }

  if (growth=="" | growth=="log") {   # log training
    if (!is.null(p.rates[["rate.log"]]))
      p.rate.log = torch$tensor(p.rates[["rate.log"]])$float()

    losses.log = x.reg$train(regr="log", p_rate=p.rate.log, steps=as.integer(steps), random_state=as.integer(random_state))
    p.log = x.reg$get_learned_params()
    ll.log = x.reg$compute_log_likelihood() %>% setNames(nm=lineages)
  }

  params = get_growth_params(timepoints_to_int,
                             rates.exp=p.exp,
                             rates.log=p.log,
                             lineages=pop_df$Lineage %>% unique(),
                             cluster=cluster)

  best = c()
  for (ll in lineages)
    if (ll.exp[ll] > ll.log[ll]) best = c(best, "exp") %>% setNames(nm=c(names(best), ll)) else
      best = c(best, "log") %>% setNames(nm=c(names(best), ll))

  return(
    params %>%
      dplyr::mutate(best_model=best[Lineage])
  )
}


get_parent_rate = function(parents, rates.df, cluster) {
  if (!cluster %in% parents$Identity) return(list("exp"=NULL, "log"=NULL))

  par = parents %>% filter(Identity==cluster) %>% dplyr::pull(Parent)  ## parent of subcl
  par.rates = rates.df %>%
    filter(Identity==par) %>%
    dplyr::select(dplyr::contains("rate.exp"), dplyr::contains("rate.log"), -dplyr::contains("p_rate")) %>%
    as.list()

  if (NA %in% unlist(par.rates)) return(list("exp"=NULL, "log"=NULL))

  return(par.rates)
}


get_growth_params = function(timepoints_to_int,
                             lineages,
                             cluster,
                             rates.exp=NULL,
                             rates.log=NULL) {
  if (!is.null(rates.exp)) params.exp = get_growth_rates_exp(rates.exp, lineages, cluster)
  if (!is.null(rates.log)) params.log = get_growth_rates_log(rates.log, lineages, cluster)

  if (is.null(rates.exp)) return(params.log)
  if (is.null(rates.log)) return(params.exp)
  return( dplyr::inner_join(params.exp, params.log, by=c("Lineage", "Identity")) )
}


get_growth_rates_exp = function(rates.df, lineages, cluster) {
  return(
    data.frame() %>%

      # initialize the columns
      tibble::add_column("fitness.exp"=NA,
                         "init_t.exp"=NA,
                         "Lineage"=as.character(NA),
                         "p_rate.exp"=NA,
                         "sigma.exp"=NA,
                         "rate.exp"=as.numeric(NA)) %>%

      # add the values
      tibble::add_row(Lineage=lineages,
                      p_rate.exp=rates.df$parent_rate,
                      fitness.exp=rates.df$fitness,
                      init_t.exp=rates.df$init_time) %>%

      # compute the rates for subclones
      dplyr::mutate(rate.exp=replace( rate.exp, !is.na(p_rate.exp), p_rate.exp * (1+fitness.exp) ),
                    Identity=cluster,
                    sigma.exp=list( setNames(object=rates.df$sigma, nm=c(0, unlist(timepoints_to_int))) )
                    ) %>%

      tibble::as_tibble()
  )
}


get_growth_rates_log = function(rates.df, lineages, cluster) {
  return(
    data.frame() %>%
      tibble::add_column("fitness.log"=NA,
                         "K.log"=NA,
                         "init_t.log"=NA,
                         "Lineage"=as.character(NA),
                         "p_rate.log"=NA,
                         "sigma.log"=NA,
                         "rate.log"=as.numeric(NA)) %>%

      tibble::add_row(fitness.log=rates.df$fitness,
                      K.log=rates.df$carr_capac,
                      init_t.log=rates.df$init_time,
                      Lineage=lineages,
                      p_rate.log=rates.df$parent_rate) %>%

      dplyr::mutate(rate.log=replace( rate.log, !is.na(p_rate.log), p_rate.log * (1+fitness.log) ),
                    Identity=cluster,
                    sigma.log=list( setNames(object=rates.df$sigma, nm=c(0, unlist(timepoints_to_int))) )
                    ) %>%

      tibble::as_tibble()
  )
}


