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
#' @return A \code{mvnmm} object with the additional tibble \code{growth.rates} containing the estimated population genetics parameters.
#'
#' @importFrom purrr is_empty
#' @importFrom dplyr select filter add_row mutate pull contains select
#'
#' @export

fit_growth_rates = function(x,
                            steps=500,
                            highlight=c(),
                            timepoints_to_int=c(),
                            growth_model="exp.log",
                            force=T,
                            tree_score=1,
                            py_pkg=NULL,
                            mutations=F) {

  highlight.cov = get_highlight(x, highlight=highlight, mutations=F)
  highlight.muts = get_highlight(x, highlight=highlight, mutations=T)
  timepoints_to_int = map_timepoints_int(x, timepoints_to_int)

  if (!mutations && have_muts_fit(x)) mutations = T

  pop_df = x %>%
    get_muller_pop(mutations=mutations, timepoints_to_int=timepoints_to_int,
                   tree_score=tree_score) %>%
    dplyr::select(-dplyr::contains("Pop.subcl"), -dplyr::contains("Pop.plot"), -dplyr::contains("theta")) %>%
    dplyr::filter(Identity %in% highlight.muts)

  rates.df = data.frame()
  evaluated = list()

  # if force is TRUE, then all the clusters are fitted again
  if (have_growth_rates(x)) {
    rates.df = x %>% get_growth_rates()
    if (force) rates.df = x %>% get_growth_rates() %>% dplyr::filter(!Identity %in% highlight.muts)

    evaluated = rates.df$Identity %>% unique()
  }

  # edges = get_muller_edges(x, mutations=mutations, tree_score=tree_score)

  for (cluster in highlight.cov) {
    if (cluster %in% evaluated) next

    cli::cli_process_start(paste0("Starting growth models inference of clone ", cluster))

    # filter the dataset
    pop_df.cl = pop_df %>%
      dplyr::filter(Identity==cluster |
                      grepl(paste0(cluster, "."), Identity, fixed=T))

    # get the identity-parents dataframe of clone "cluster"
    # parents = get_parents(x, highlight=cluster, tree_score=tree_score)
    parents = pop_df.cl %>%
      dplyr::select(Parent, Identity) %>% unique()

    rates.df = fit_growth_utils(rates.df=rates.df,
                                cluster=cluster,
                                pop_df.cl=pop_df.cl,
                                parents=parents,
                                growth_model=growth_model,
                                steps=steps,
                                timepoints_to_int=timepoints_to_int,
                                py_pkg=py_pkg)

    evaluated = c(cluster, evaluated)

    cli::cli_process_done()
  }

  x = add_growth_rates(x, rates.df)
  x = add_tp_int(x, timepoints_to_int)

  return(x)
}


fit_growth_utils = function(rates.df,
                            cluster,
                            pop_df.cl,
                            parents,
                            growth_model,
                            steps,
                            timepoints_to_int,
                            py_pkg) {

  # first in the clonal cluster of ISs
  rates.df = rates.df %>%
    fit_growth_clones(clusters=cluster,
                      pop_df.cl=pop_df.cl,
                      parents=parents,
                      growth_model=growth_model,
                      steps=steps,
                      clonal=TRUE,
                      timepoints_to_int=timepoints_to_int,
                      py_pkg=py_pkg)

  subclones = c(parents$Identity, parents$Parent) %>%
    sort_clusters_edges(edges=parents) # get the clusters sorted by tree hierarchy

  # fit the growth rate in all the parents
  rates.df = rates.df %>%
    fit_growth_clones(clusters=subclones,
                      pop_df.cl=pop_df.cl,
                      parents=parents,
                      growth_model=growth_model,
                      steps=steps,
                      clonal=FALSE,
                      timepoints_to_int=timepoints_to_int,
                      py_pkg=py_pkg)

  return(rates.df)
}


fit_growth_clones = function(rates.df,
                             clusters,
                             pop_df.cl,
                             parents,
                             growth_model,
                             steps,
                             clonal,
                             timepoints_to_int,
                             py_pkg=NULL) {

  if (!clonal & nrow(parents)==0) return(rates.df)

  # if (!clonal) clusters = sort_clusters_edges(clusters, parents)

  for (subcl in clusters)
    # par = pop_df.cl %>% dplyr::filter(Identity==subcl) %>% dplyr::pull(Parent) %>% unique()
    if ( (!subcl %in% rates.df$Identity) && (subcl %in% pop_df.cl$Identity) && (subcl != "P") )
        rates.df = run_py_growth(rates.df,
                                 pop_df.cl=pop_df.cl,
                                 cluster=subcl,
                                 timepoints_to_int=timepoints_to_int,
                                 # p.rates=get_parent_rate(par, rates.df, subcl),
                                 clonal=clonal,
                                 growth_model=growth_model,
                                 steps=steps,
                                 py_pkg=py_pkg)

  return(rates.df)
}


run_py_growth = function(rates.df,
                         pop_df.cl,
                         cluster,
                         timepoints_to_int,
                         steps=500,
                         # p.rates=list("exp"=NULL, "log"=NULL),
                         clonal=FALSE,
                         growth_model="exp.log",
                         random_state=25,
                         py_pkg=NULL) {

  par = pop_df.cl %>% dplyr::filter(Identity==cluster) %>% dplyr::pull(Parent) %>% unique()
  p.rates = get_parent_rate(par, rates.df, cluster)

  torch = reticulate::import("torch")
  if (is.null(py_pkg)) py_pkg = reticulate::import("pylineaGT")

  input.clone = pop_df.cl %>%
    dplyr::filter(Identity == cluster) %>%

    ## Run with population size
    # dplyr::select(-Frequency) %>%
    # dplyr::mutate(Population=ifelse((Generation == 0 & clonal), 1, Population)) %>%
    # tidyr::pivot_wider(id_cols="Generation", names_from="Lineage", values_from="Population") %>%

    ## Run with re-scaled population size
    dplyr::mutate(Frequency=Frequency * 1000) %>%
    dplyr::select(-Population) %>%
    dplyr::mutate(Frequency=ifelse((Generation == 0 & clonal), 1, Frequency)) %>%
    tidyr::pivot_wider(id_cols="Generation", names_from="Lineage", values_from="Frequency") %>%

    dplyr::arrange(Generation)

  times = input.clone$Generation
  y = input.clone[2:ncol(input.clone)]
  lineages = colnames(y)

  times = torch$tensor(times)$int()$unsqueeze(as.integer(1))
  y = torch$tensor(y %>% as.matrix())

  # Remove ######
  print(times)
  print(y)

  p.rate.exp = p.rate.log = NULL

  x.reg = py_pkg$explogreg$Regression(times, y)
  if (grepl("exp", growth_model)) {  # exp training
    if (!is.null(p.rates[["exp"]])) p.rate.exp = torch$tensor(p.rates[["exp"]])$float()

    losses.exp = x.reg$train(regr="exp", p_rate=p.rate.exp, steps=as.integer(steps), random_state=as.integer(random_state))
    p.exp = x.reg$get_learned_params()
    ll.exp = x.reg$compute_log_likelihood() %>% setNames(nm=lineages)

    # Remove ######
    print(p.exp)
  }

  if (grepl("log", growth_model)) {   # log training
    if (!is.null(p.rates[["log"]])) p.rate.log = torch$tensor(p.rates[["log"]])$float()

    losses.log = x.reg$train(regr="log", p_rate=p.rate.log, steps=as.integer(steps), random_state=as.integer(random_state))
    p.log = x.reg$get_learned_params()
    ll.log = x.reg$compute_log_likelihood() %>% setNames(nm=lineages)
  }

  params = get_growth_params(timepoints_to_int,
                             rates.exp=p.exp,
                             rates.log=p.log,
                             lineages=pop_df.cl$Lineage %>% unique(),
                             cluster=cluster)

  best = c()
  for (ll in lineages)
    if (ll.exp[ll] > ll.log[ll])
      best = c(best, "exp") %>% setNames(nm=c(names(best), ll)) else
      best = c(best, "log") %>% setNames(nm=c(names(best), ll))


  if (purrr::is_empty(rates.df))
    return(
      params %>%
        dplyr::mutate(best_model=best[Lineage])
    )

  return(
    rates.df %>%
      dplyr::add_row(
        params %>%
          dplyr::mutate(best_model=best[Lineage])
      )
  )
}

