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
                            tree_score=1,
                            py_pkg=NULL) {

  highlight.cov = get_highlight(x, highlight=highlight, mutations=F)
  highlight.muts = get_highlight(x, highlight=highlight, mutations=T)
  timepoints_to_int = map_timepoints_int(x, timepoints_to_int)

  if (have_muts_fit(x)) mutations = T else mutations = F

  pop_df = x %>%
    get_muller_pop(mutations=mutations, timepoints_to_int=timepoints_to_int,
                   tree_score=tree_score)

  rates.df = data.frame()
  evaluated = list()

  # if force is TRUE, then all the clusters are fitted again
  if (force & have_growth_rates(x)) {
    rates.df = x %>% get_growth_rates() %>% dplyr::filter(!Identity %in% highlight.muts)
    evaluated = rates.df$Identity %>% unique()
  } else if (!force & have_growth_rates(x)) {
    rates.df = x %>% get_growth_rates()
    evaluated = rates.df$Identity %>% unique()
  }

  edges = get_muller_edges(x, mutations=mutations, tree_score=tree_score)

  for (cluster in highlight.cov) {
    if (!cluster %in% evaluated) {

      cli::cli_process_start(paste0("Starting growth models inference of clone ", cluster))

      # filter the dataset
      pop_df.cl = pop_df %>%
        dplyr::filter(grepl(paste(cluster,".",sep=""), Identity) | Identity == cluster)

      # get the identity-parents dataframe of clone "cluster"
      # parents = get_parents(x, highlight=cluster, tree_score=tree_score)
      parents = edges %>%
        dplyr::filter(grepl(paste(cluster,".",sep=""), Parent) |
                        grepl(paste(cluster,".",sep=""), Identity) |
                        Identity==cluster |
                        Parent==cluster)

      # first in the clonal cluster of ISs
      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=cluster,
                                   pop_df=pop_df.cl,
                                   parents=parents,
                                   growth=growth_model,
                                   steps=steps,
                                   clonal=TRUE,
                                   timepoints_to_int=timepoints_to_int,
                                   py_pkg=py_pkg)

      # fit the growth rate in all the parents
      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=parents$Parent %>% unique(),
                                   pop_df=pop_df.cl,
                                   parents=parents,
                                   growth=growth_model,
                                   steps=steps,
                                   clonal=F,
                                   timepoints_to_int=timepoints_to_int,
                                   py_pkg=py_pkg)

      # then in all the children
      rates.df = rates.df %>%
        fit_growth_multiple_clones(clusters=parents$Identity %>% unique(),
                                   pop_df=pop_df.cl,
                                   parents=parents,
                                   growth=growth_model,
                                   steps=steps,
                                   clonal=F,
                                   timepoints_to_int=timepoints_to_int,
                                   py_pkg=py_pkg)

      cli::cli_process_done()

      evaluated = c(cluster, evaluated)
    }
  }

  x = add_growth_rates(x, rates.df)
  x = add_tp_int(x, timepoints_to_int)

  return(x)
}


fit_growth_multiple_clones = function(rates.df,
                                      clusters,
                                      pop_df,
                                      parents,
                                      growth,
                                      steps,
                                      clonal,
                                      timepoints_to_int,
                                      py_pkg=NULL) {

  if (!clonal & nrow(parents)==0) return(rates.df)

  # get the clusters sorted by tree hierarchy
  if (!clonal) clusters = sort_clusters_edges(clusters, parents)

  for (subcl in clusters) {

    if (!subcl %in% (rates.df$Identity) && (subcl %in% pop_df$Identity) && (subcl != "P")) {
      if (purrr::is_empty(rates.df))
        rates.df = fit_growth_single_clone(pop_df=pop_df,
                                           cluster=subcl,
                                           timepoints_to_int=timepoints_to_int,
                                           p.rates=get_parent_rate(parents, rates.df, subcl),
                                           clonal=clonal,
                                           growth=growth,
                                           steps=steps,
                                           py_pkg=py_pkg)
      else rates.df = rates.df %>%
          dplyr::add_row(
            fit_growth_single_clone(pop_df=pop_df,
                                    cluster=subcl,
                                    timepoints_to_int=timepoints_to_int,
                                    p.rates=get_parent_rate(parents, rates.df, subcl),
                                    clonal=clonal,
                                    growth=growth,
                                    steps=steps,
                                    py_pkg=py_pkg)
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
    if (!is.null(p.rates[["exp"]]))
      p.rate.exp = torch$tensor(p.rates[["exp"]])$float()

    losses.exp = x.reg$train(regr="exp", p_rate=p.rate.exp, steps=as.integer(steps), random_state=as.integer(random_state))
    p.exp = x.reg$get_learned_params()
    ll.exp = x.reg$compute_log_likelihood() %>% setNames(nm=lineages)
  }

  if (growth=="" | growth=="log") {   # log training
    if (!is.null(p.rates[["log"]]))
      p.rate.log = torch$tensor(p.rates[["log"]])$float()

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

