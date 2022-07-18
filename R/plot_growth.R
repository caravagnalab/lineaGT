#' Exponential fitting
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
#' @export plot_exp_fit

plot_exp_fit = function(x,
                        highlight=c(),
                        min_frac=0,
                        mutations=T,
                        label="",
                        timepoints_to_int=list(),
                        fit=F) {

  if (purrr::is_empty(timepoints_to_int)) timepoints_to_int = map_timepoints_int(x, timepoints_to_int)
  highlight = get_highlight(x, min_frac, highlight, mutations=mutations)
  color_palette = highlight_palette(x, highlight, label)

  if (fit)
    x = fit_growth(x, highlight=highlight, timepoints_to_int=timepoints_to_int, force=F)

  highlight = intersect(highlight, x %>% get_growth_rates() %>% dplyr::pull(Identity))

  pop_df = get_muller_pop(x, mutations=mutations, label=label, timepoints_to_int=timepoints_to_int) %>%
    filter(Identity %in% highlight)
  regr.df = get_regression_df(x, highlight=highlight)

  p = pop_df %>%
    ggplot() +
    geom_point(aes(x=Generation, y=Population, color=Lineage)) +
    geom_line(data=regr.df, aes(x=x, y=y.exp, color=Lineage), linetype="dashed") +
    geom_line(data=regr.df, aes(x=x, y=y.log, color=Lineage), linetype="solid") +
    geom_vline(data=regr.df, aes(xintercept=init_t.exp), linetype="dashed", size=.1) +
    geom_vline(data=regr.df, aes(xintercept=init_t.log), linetype="solid", size=.1) +
    facet_grid(rows=vars(Identity), cols=vars(Lineage), scales="free_y") +
    my_ggplot_theme()

  return(p)
}


#' Exponential rate
#'
#' @description Function to visualize the growth coefficients for each clone and lineage.
#'
#' @param x a mvnmm object.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param mutations Boolean. If set to \code{TRUE}, the growth will be visualize for each cluster of mutations.
#' @param label a character corresponding to the label of the run to visualize.
#' @param timepoints_to_int a list to map each \code{timepoint} value to an integer.
#'
#' @examples
#' if (FALSE) plot_exp_rate(x)
#'
#' @import ggplot2
#'
#' @export plot_exp_rate

plot_exp_rate = function(x,
                         highlight=c(),
                         min_frac=0,
                         mutations=F,
                         label="",
                         timepoints_to_int=list()) {

  if (purrr::is_empty(timepoints_to_int)) timepoints_to_int = map_timepoints_int(x)

  highlight = get_highlight(x, min_frac, highlight, mutations=mutations, label=label) %>%
    color_palette = highlight_palette(x, highlight, label)

  highlight = highlight %>% stringr::str_replace("C_||C", "")
  names(color_palette) = names(color_palette) %>% stringr::str_replace("C_||C", "")

  p = x %>%
    get_muller_pop(timepoints_to_int=timepoints_to_int, mutations=mutations, label=label) %>%
    dplyr::mutate(Identity=Identity %>% stringr::str_replace("C_||C","")) %>%
    dplyr::filter(Identity %in% highlight) %>%
    dplyr::arrange(lm_r) %>%
    ggplot(aes(x=Identity, y=lm_r, ymax=lm_r, ymin=0, color=Identity)) +
    geom_linerange(alpha=.3) +
    geom_point(alpha=.3) +
    facet_wrap(~Lineage) +
    scale_color_manual(values=color_palette[highlight]) +
    my_ggplot_theme() +
    theme(aspect.ratio=1) +
    xlab("Clusters") + ylab("Exp rate") + labs(color="Clusters")
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
    dplyr::select(-dplyr::contains("rate"), -dplyr::contains("K.log")) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(x, Identity)

  return(regr_df)
}
