#' Visualize the regression given the infered growth rates.
#'
#' @description Function to visualize the growht of each clone and lineage.
#'
#' @param x a mvnmm object.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param mutations Boolean. If set to \code{TRUE}, the growth will be visualize for each cluster of mutations.
#' @param label a character corresponding to the label of the run to visualize.
#' @param timepoints_to_int a list to map each \code{timepoint} value to an integer.
#' @param fit add
#'
#' @examples
#' if (FALSE) plot_exp_fit(x)
#'
#' @import ggplot2
#'
#' @export plot_growth_regression

plot_growth_regression = function(x,
                        highlight=c(),
                        min_frac=0,
                        mutations=F,
                        label="",
                        timepoints_to_int=list(),
                        fit=F) {

  if (purrr::is_empty(timepoints_to_int)) timepoints_to_int = map_timepoints_int(x, timepoints_to_int)
  highlight = get_highlight(x, min_frac, highlight, mutations=mutations)
  color_palette = highlight_palette(x, highlight, label)

  if (fit)
    x = fit_growth_rates(x, highlight=highlight, timepoints_to_int=timepoints_to_int, force=fit)

  # keep only the clusters with growth model fitted
  highlight = intersect(highlight, x %>% get_growth_rates() %>% dplyr::pull(Identity))

  pop_df = get_muller_pop(x, mutations=mutations, label=label, timepoints_to_int=timepoints_to_int) %>%
    filter(Identity %in% highlight)
  regr.df = get_regression_df(x, highlight=highlight) %>%
    tidyr::pivot_longer(starts_with("y."), names_to=c("else","type"), names_sep="[.]", values_to="y", names_repair="minimal") %>%
    dplyr::mutate(sigma=ifelse(type=="log", sigma.log, sigma.exp),
                  init_t=ifelse(type=="log", init_t.log, init_t.exp),
                  rate=ifelse(type=="log", rate.log, rate.exp),
                  K=ifelse(type=="log", K.log, NA)) %>%
    dplyr::mutate(y.min=ifelse(type=="log",
                               K / ( 1 + (K-1) * exp( (-rate*(x-init_t)) - sigma ) ),
                               exp( rate * (x-init_t) - sigma ))) %>%
    dplyr::mutate(y.max=ifelse(type=="log",
                               K / ( 1 + (K-1) * exp( (-rate*(x-init_t)) + sigma ) ),
                               exp( rate * (x-init_t) + sigma ))) %>%
    dplyr::select(-"else",-dplyr::ends_with("exp"),-dplyr::ends_with("log")) %>%
    dplyr::mutate(type=ifelse(type=="log", "Logistic", "Exponential"))


  color_palette = c("firebrick","steelblue"); names(color_palette) = c("Exponential","Logistic")
  color_palette = c(color_palette, x$color_palette)
  p = pop_df %>%
    ggplot() +
    geom_point(aes(x=Generation, y=Population), alpha=.5, size=.7) +
    geom_line(data=regr.df, aes(x=x, y=y, color=type), size=.6, alpha=.8) +
    geom_vline(data=regr.df, aes(xintercept=init_t, color=type), linetype="dashed", size=.2, alpha=.5) +
    geom_errorbar(data=regr.df, aes(x=x, y=y, ymin=y.min, ymax=y.max, color=type), width=.5, position=position_dodge(width=0.5)) +
    facet_grid(rows=vars(Identity), cols=vars(Lineage), scales="free_y") +
    scale_color_manual(values=color_palette, breaks=c("Exponential","Logistic")) +
    labs(color="") +
    my_ggplot_theme()

  return(p)
}


#' Visualize the infered growth rates.
#'
#' @description Function to visualize the growth coefficients for each clone and lineage.
#'
#' @param x a mvnmm object.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param mutations Boolean. If set to \code{TRUE}, the growth will be visualize for each cluster of mutations.
#' @param label a character corresponding to the label of the run to visualize.
#' @param timepoints_to_int a list to map each \code{timepoint} value to an integer.
#' @param fit add
#'
#' @examples
#' if (FALSE) plot_growth_rates(x)
#'
#' @import ggplot2
#'
#' @export plot_growth_rates

plot_growth_rates = function(x,
                         highlight=c(),
                         min_frac=0,
                         mutations=F,
                         label="",
                         timepoints_to_int=list(),
                         fit=F) {

  if (purrr::is_empty(timepoints_to_int)) timepoints_to_int = map_timepoints_int(x, timepoints_to_int=timepoints_to_int)
  highlight = get_highlight(x, min_frac, highlight, mutations=mutations)

  if (fit)
    x = fit_growth_rates(x, highlight=highlight, timepoints_to_int=timepoints_to_int, force=F)

  # keep only the clusters with growth model fitted
  highlight = intersect(highlight, x %>% get_growth_rates() %>% dplyr::pull(Identity))

  rates = x %>% get_growth_rates() %>% dplyr::select(dplyr::starts_with("rate"), Identity, Lineage) %>%
    tidyr::pivot_longer(cols=dplyr::starts_with("rate"), names_to="type", names_prefix="rate.", values_to="rate") %>%
    dplyr::mutate(type=ifelse(type=="exp", "Exponential", "Logistic"))

  color_palette = c("firebrick","steelblue"); names(color_palette) = c("Exponential","Logistic")
  p = rates %>%
    # dplyr::select(-dplyr::contains(select))
    # dplyr::mutate(Identity=Identity %>% stringr::str_replace("C_||C","")) %>%
    dplyr::filter(Identity %in% highlight) %>%
    # dplyr::arrange(lm_r) %>%
    ggplot(aes(x=Identity, y=rate, ymax=rate, ymin=0, color=type)) +
    geom_linerange(alpha=.8, position=position_dodge2(width=.5)) +
    geom_point(alpha=.8, position=position_dodge2(width=.5)) +
    facet_wrap(~Lineage) +
    scale_color_manual(values=color_palette) +
    my_ggplot_theme() +
    # theme(aspect.ratio=1) +
    xlab("Clusters") + ylab("Growth rate") + labs(color="")

  return(p)
}




# function to construct a regression dataframe per clone/lineage
get_regression_df = function(x, highlight) {
  pop_df = x %>% get_muller_pop(mutations=TRUE, exp_coef=F, timepoints_to_int=get_tp_to_int(x)) %>%
    filter(Identity %in% highlight)

  rates = x %>% get_growth_rates() %>%
    dplyr::filter(Identity %in% highlight) %>%
    dplyr::select(dplyr::contains("rate"),
                  dplyr::contains("K.log"),
                  dplyr::contains("init_t"),
                  dplyr::contains("sigma"),
                  "Lineage", "Identity")

  tmax = max(pop_df$Generation)
  tmin = 0

  n_lins = rates$Lineage %>% unique() %>% length()
  n_cls = length(highlight)

  regr_df = data.frame(x=rep( tmin:tmax, times=n_lins*n_cls ),
                       Lineage=rep( rates$Lineage %>% unique(), each=(tmax-tmin+1)*n_cls ),
                       Identity=rep( rates$Identity %>% unique(), times=(tmax-tmin+1)*n_lins) ) %>%
    dplyr::inner_join(rates, by=c("Lineage", "Identity")) %>%
    dplyr::mutate(y.exp=exp( rate.exp * (x-init_t.exp) ),
                  y.log=K.log / ( 1 + (K.log-1) * exp( -rate.log*(x-init_t.log) ) )) %>%
    dplyr::select(-dplyr::contains("p_rate")) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(x, Identity) %>%
    dplyr::mutate(sigma.exp=unlist(sigma.exp)[as.character(x)], sigma.log=unlist(sigma.log)[as.character(x)])


  return(regr_df)
}
