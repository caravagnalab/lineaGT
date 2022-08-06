#' Visualize the regression given the infered growth rates.
#'
#' @description Function to visualize the growht of each clone and lineage.
#'
#' @param x a mvnmm object.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param mutations Boolean. If set to \code{TRUE}, the growth will be visualize for each cluster of mutations.
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
                                  timepoints_to_int=list(),
                                  fit=F) {

  timepoints_to_int = map_timepoints_int(x, timepoints_to_int)
  highlight = get_highlight(x, min_frac, highlight, mutations=mutations)
  color_palette = highlight_palette(x, highlight)

  if (fit)
    x = fit_growth_rates(x, highlight=highlight, timepoints_to_int=timepoints_to_int, force=fit)
  else if (purrr::is_empty(get_growth_rates(x)))
    x = fit_growth_rates(x, highlight=highlight, timepoints_to_int=timepoints_to_int, force=T)

  # keep only the clusters with growth model fitted
  highlight = intersect(highlight, x %>% get_growth_rates() %>% dplyr::pull(Identity))

  pop_df = x %>%
    get_muller_pop(mutations=mutations, timepoints_to_int=timepoints_to_int) %>%
    filter(Identity %in% highlight)

  regr.df = get_regression_df(x, pop_df=pop_df, highlight=highlight)

  color_palette = c("firebrick","steelblue"); names(color_palette) = c("Exponential","Logistic")
  color_palette = c(color_palette, x$color_palette)
  return(
    pop_df %>%
    ggplot() +
    geom_point(aes(x=Generation, y=Population), alpha=.5, size=.7) +
    geom_line(data=regr.df, aes(x=x, y=y, color=type), size=.6, alpha=.8) +
    geom_vline(data=regr.df, aes(xintercept=init_t, color=type), linetype="dashed", size=.2, alpha=.5) +
    geom_errorbar(data=regr.df, aes(x=x, y=y, ymin=y.min, ymax=y.max, color=type), width=.5, position=position_dodge(width=0.5)) +
    facet_grid(rows=vars(Identity), cols=vars(Lineage), scales="free_y") +
    scale_color_manual(values=color_palette, breaks=c("Exponential","Logistic")) +
    labs(color="") + my_ggplot_theme()
  )
  return(p)
}


# function to construct a regression dataframe per clone/lineage
get_regression_df = function(x, pop_df, highlight) {
  rates = x %>% get_growth_rates() %>%
    dplyr::filter(Identity %in% highlight)
  tmax = max(pop_df$Generation)
  tmin = 0

  n_lins = rates$Lineage %>% unique() %>% length()
  n_cls = length(highlight)

  regr_df = rates %>%
    dplyr::select(-dplyr::starts_with("p_rate"), -dplyr::starts_with("fitness")) %>%
    dplyr::mutate(x=list(tmin:tmax)) %>%
    tidyr::unnest(x) %>%
    dplyr::mutate(sigma.exp=unlist(sigma.exp)[as.character(x)],
                  sigma.log=unlist(sigma.log)[as.character(x)]) %>%
    dplyr::arrange(x, Identity) %>%
    # dplyr::select(-dplyr::starts_with("rate"), -dplyr::starts_with("K")) %>%

    tidyr::pivot_longer(starts_with("rate."), names_to=c("else","type"), names_sep="[.]",
                        values_to="rate", names_repair="minimal") %>%

    dplyr::mutate(init_t=ifelse(type=="log", init_t.log, init_t.exp),
                  args=ifelse(type=="log", -rate*(x-init_t), rate*(x-init_t) ),
                  y=ifelse(type=="log", K.log/( 1+(K.log-1)*exp(args) ), exp(args)),
                  y.min=ifelse(type=="log", K.log/( 1+(K.log-1)*exp(args-sigma.log) ), exp(args-sigma.exp)),
                  y.max=ifelse(type=="log", K.log/( 1+(K.log-1)*exp(args+sigma.log) ), exp(args+sigma.exp)) ) %>%
    dplyr::select(Lineage, Identity, init_t, type, x, y, y.min, y.max) %>%
    dplyr::mutate(type=ifelse(type=="log", "Logistic", "Exponential"))

  return(regr_df)
}



#' Visualize the infered growth rates.
#'
#' @description Function to visualize the growth coefficients for each clone and lineage.
#'
#' @param x a mvnmm object.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param mutations Boolean. If set to \code{TRUE}, the growth will be visualize for each cluster of mutations.
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
