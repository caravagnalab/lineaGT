compute_rates = function(x, cluster, timepoints_to_int=list()) {
  pop = x %>%
    get_muller_pop(mutations=T, timepoints_to_int=timepoints_to_int) %>%
    filter(grepl(cluster, Identity)) %>%
    dplyr::select(-starts_with("lm"))
  parents = get_parents(x, highlight=cluster)

  rates = list(); regr.df = data.frame()
  rates[[cluster]] = run_exp_single_clone(pop, cluster, clonal=TRUE)
  regr.df = rbind(regr.df, get_regression_df(pop, rates[[cluster]], cluster))

  evaluated = c(cluster)
  for (pp in (parents$Parent %>% unique())) {
    ## check if "pp" has a parent different from "cluster"
    if (pp != cluster && !pp %in% evaluated) {
      print(pp)
      pp.par = parents %>% filter(Identity==pp) %>% dplyr::pull(Parent)  ## parent of pp

      parent.rate = compute_parent_rate(rates[[pp.par]])
      ## case2 -> it is another subclone
      if (pp.par != cluster) {
        ## TODO add check if we have not computed yet its parent
        print("")
      }

      print(parent.rate$exp$rate)
      rates[[pp]] = run_exp_single_clone(pop, pp, clonal=FALSE, p.rates=parent.rate$exp$rate)
      regr.df = rbind(regr.df, get_regression_df(pop, rates[[pp]], pp))

      evaluated = c(evaluated, pp)
    }
  }

  for (pp in (parents$Identity %>% unique())) {
    ## check if pp has not been evaluated yet
    if (!pp %in% evaluated) {
      pp.par = parents %>% filter(Identity==pp) %>% dplyr::pull(Parent)  ## parent of pp
      parent.rate = compute_parent_rate(rates[[pp.par]])

      rates[[pp]] = run_exp_single_clone(pop, pp, clonal=FALSE, p.rates=parent.rate$exp$rate)
      regr.df = rbind(regr.df, get_regression_df(pop, rates[[pp]], pp))
      evaluated = c(evaluated, pp)
    }
  }

  return(list("rates"=rates, "df"=regr.df))
}


run_exp_single_clone = function(pop_df, cluster, p.rates=list(), clonal=FALSE) {
  if (clonal) {
    input.clone = pop_df %>%
      filter(Identity == cluster) %>%
      dplyr::select(-Frequency) %>%
      ## mutate first value to 1 only for the CLONAL cluster
      dplyr::mutate(Population=ifelse(Generation == 0, 1, Population)) %>%
      pivot_wider(id_cols="Generation", names_from="Lineage", values_from="Population")
  } else {
    input.clone = pop_df %>%
      filter(Identity == cluster) %>%
      dplyr::select(-Frequency) %>%
      pivot_wider(id_cols="Generation", names_from="Lineage", values_from="Population")
  }


  time = input.clone$Generation; time = torch$tensor(time)$int()$unsqueeze(as.integer(1))
  y = input.clone[2:ncol(input.clone)] %>% as.matrix(); y = torch$tensor(y)

  parent_rate = NULL
  if (!purrr::is_empty(p.rates)) parent_rate = torch$tensor(p.rates)$float()

  print(parent_rate)
  print(cluster)

  x.reg = Regression(time, y, p_rate=parent_rate)

  ## exp training
  losses.exp = x.reg$train(regr="exp")
  rr.exp = x.reg$get_learned_params()
  rates.exp = rr.exp$AutoDelta.fitness$detach()$numpy(); names(rates.exp) = pop$Lineage %>% unique()

  ## log training
  # x.reg$p_rate = NULL
  # losses.log = x.reg$train(regr="log")
  # rr.log = x.reg$get_learned_params()
  # rates.log = rr.log$AutoDelta.rate$detach()$numpy(); names(rates.log) = pop$Lineage %>% unique()
  # carr_cap.log = rr.log$AutoDelta.carr_capac$detach()$numpy(); names(carr_cap.log) = pop$Lineage %>% unique()

  return(list(
      exp=list("fitness"=rates.exp, "p.rate"=p.rates %>% unlist())
      # "log"=list("fitness"=rates.log, "p.rate"=parent_rate, "carrying_cap"=carr_cap.log)
  ))
}


compute_parent_rate = function(rates) {
  if (is.null(rates$exp[["p.rate"]])) rates.exp = rates$exp[["fitness"]]
  else rates.exp = rates$exp[["p.rate"]] * (1 + rates$exp[["fitness"]])

  # if (is.null(rates["log"]["p.rate"])) rates.log = rates["log"]["fitness"]
  # else rates.log = rates["log"]["p.rate"] * (1 + rates["log"]["fitness"])
  # carr_cap.log = rates["log"]["carrying_cap"]

  return(
    list(
      exp=list("rate"=rates.exp)
      # log=list("rate"=rates.log, "carr_cap"=carr_cap.log)
      )
    )
}


get_regression_df = function(pop_df, rates, cluster) {
  # if (is.null(rates["exp"]["p.rate"])) rates.exp = rates["exp"]["fitness"]
  # else rates.exp = rates["exp"]["p.rate"] * (1 + rates["exp"]["fitness"])

  if (is.null(rates$exp[["p.rate"]])) rates.exp = rates$exp[["fitness"]]
  else rates.exp = rates$exp[["p.rate"]] * (1 + rates$exp[["fitness"]])

  # if (is.null(rates["log"]["p.rate"])) rates.log = rates["log"]["fitness"]
  # else rates.log = rates["log"]["p.rate"] * (1 + rates["log"]["fitness"])
  # carr_cap.log = rates["log"]["carrying_cap"]

  regr_df = data.frame(x=rep(1:280, times=(pop_df$Lineage %>% unique() %>% length())),
                       Lineage=rep(pop_df$Lineage %>% unique(), each=280)) %>%
    mutate(y.exp=exp(rates.exp[Lineage] * x)) %>%
    # mutate(y.log=carr_cap.log[Lineage] / (1 + (carr_cap.log[Lineage]) * exp(-rates.log[Lineage] * x))) %>%
    mutate(Identity=cluster) %>%
    tibble::as_tibble()

  return(regr_df)
}
