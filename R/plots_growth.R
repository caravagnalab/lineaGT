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
                                  fit=F,
                                  show_best=T,
                                  ratio=NULL) {

  timepoints_to_int = map_timepoints_int(x, timepoints_to_int)
  highlight = get_highlight(x, min_frac=min_frac, highlight=highlight, mutations=mutations)
  color_palette = highlight_palette(x, highlight)

  if (fit)
    x = fit_growth_rates(x, highlight=highlight, timepoints_to_int=timepoints_to_int, force=fit)
  else if (purrr::is_empty(get_growth_rates(x)))
    x = fit_growth_rates(x, highlight=highlight, timepoints_to_int=timepoints_to_int, force=T)

  # keep only the clusters with growth model fitted
  highlight = intersect(highlight, x %>% get_growth_rates() %>% dplyr::pull(Identity))

  pop_df = x %>%
    get_muller_pop(mutations=mutations, timepoints_to_int=timepoints_to_int) %>%
    dplyr::filter(Identity %in% highlight)

  regr.df = get_regression_df(x, pop_df=pop_df, highlight=highlight)

  color_palette = c("firebrick","steelblue"); names(color_palette) = c("Exponential","Logistic")
  color_palette = c(color_palette, x$color_palette)

  pl = pop_df %>%

    ggplot() +

    # geom_point(aes(x=Generation, y=Population), alpha=.5, size=.7) +
    geom_point(aes(x=Generation, y=Frequency), alpha=.5, size=.7) +

    geom_line(data=filter(regr.df, type==best_model), aes(x=x, y=y/1000, color=type), size=.7, alpha=.9) +

    geom_vline(data=filter(regr.df, type==best_model),
               aes(xintercept=init_t, color=type), linetype="dashed", linewidth=.3, alpha=.7) +

    geom_errorbar(data=filter(regr.df, type==best_model), 
                  aes(x=x, y=y/1000, ymin=y.min/1000, ymax=y.max/1000, color=type),
                  width=.5, position=position_dodge(width=0.5), alpha=.9, size=.6) +
    facet_grid(rows=vars(Identity), cols=vars(Lineage), scales="free_y") +
    scale_color_manual(values=color_palette, breaks=unique(filter(regr.df, type==best_model)$best_model)) +
    labs(color="") + my_ggplot_theme()

  if (!is.null(ratio))
    pl = pl + theme(aspect.ratio=ratio)

  if (!show_best)
    return(
      pl +
        geom_line(data=filter(regr.df, type!=best_model), aes(x=x, y=y/1000, color="gainsboro"), size=.4, alpha=.7) +
        geom_vline(data=filter(regr.df, type!=best_model),
                   aes(xintercept=init_t, color="gainsboro"), linetype="dashed", size=.1, alpha=.4) +
        geom_errorbar(data=filter(regr.df, type!=best_model), 
                      aes(x=x, y=y/1000, ymin=y.min/1000, ymax=y.max/1000, color="gainsboro"),
                      width=.5, position=position_dodge(width=0.5), alpha=.7, size=.4) +
        scale_color_manual(values=color_palette, breaks=unique(regr.df$best_model))
    )

  return( pl )
}



# function to construct a regression dataframe per clone/lineage
get_regression_df = function(x, pop_df, highlight) {
  rates = x %>% get_growth_rates() %>%
    dplyr::filter(Identity %in% highlight)

  tmax = max(pop_df$Generation)
  tmin = 0

  n_lins = rates$Lineage %>% unique() %>% length()
  n_cls = length(highlight)

  compute_credint = function(post_samples, p_rate, x_values, type) {
    if (is.na(post_samples)) return(NA)
    post_fitness = post_samples$post_fitness[[1]]
    post_rate = dplyr::case_when(
      is.null(p_rate) | is.na(p_rate) ~ post_fitness,
      .default = p_rate * (1+post_fitness)
    )
    post_init_t = dplyr::case_when(
      is.na(post_samples$post_init_time[[1]]) ~ 0,
      .default = post_samples$post_init_time[[1]]
    )

    if (type == "exp") {
      args = post_rate * (x_values - post_init_t)
      return( exp( args ) )
    }

    post_carr_capacity = post_samples$post_carr_capacity[[1]]
    args = - post_rate * (x_values - post_init_t)
    return( post_carr_capacity / ( 1 + (post_carr_capacity - 1) * exp(args) ) )
  }

  if (!"posterior_samples" %in% colnames(rates)) rates = rates %>% dplyr::mutate(posterior_samples=NA)

  regr_df = rates %>%
    dplyr::mutate(x=list(tmin:tmax)) %>%
    tidyr::unnest(x) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      args=dplyr::case_when(
        type=="log" ~ -rate*(x-init_t),
        type=="exp" ~ rate*(x-init_t)
        ),
      y=dplyr::case_when(
        type=="log" ~ K / ( 1 + (K-1)*exp(args) ),
        type=="exp" ~ exp(args)
      )
    ) %>%
    dplyr::mutate(y_credint=list(compute_credint(posterior_samples, p_rate, x, type))) %>%

    dplyr::mutate(y.min=replace(NA, as.character(x) %in% names(sigma), y-sigma[as.character(x)]),
                  y.max=replace(NA, as.character(x) %in% names(sigma), y+sigma[as.character(x)])) %>%

    dplyr::mutate(y.min=replace(y.min, !is.na(y_credint), quantile(y_credint, 0.05, na.rm=TRUE)),
                  y.max=replace(y.max, !is.na(y_credint), quantile(y_credint, 0.95, na.rm=TRUE)),
                  y_credint=replace(y_credint, is.na(y_credint), y)) %>%

    dplyr::ungroup() %>%

    dplyr::select(Lineage, Identity, type, best_model, init_t, x, y, y.min, y.max) %>%
    dplyr::mutate(type=ifelse(type=="log", "Logistic", "Exponential")) %>%
    dplyr::mutate(best_model=ifelse(best_model=="log", "Logistic", "Exponential")) %>%
    dplyr::mutate(Identity=factor(Identity, levels=highlight))

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
                         fit=F,
                         show_best=T) {

  if (purrr::is_empty(timepoints_to_int)) timepoints_to_int = map_timepoints_int(x, timepoints_to_int=timepoints_to_int)
  highlight = get_highlight(x, min_frac=min_frac, highlight=highlight, mutations=mutations)

  if (fit)
    x = fit_growth_rates(x, highlight=highlight, timepoints_to_int=timepoints_to_int, force=F)

  # keep only the clusters with growth model fitted
  highlight = intersect(highlight, x %>% get_growth_rates() %>% dplyr::pull(Identity))

  rates = x %>%
    get_growth_rates() %>%
    dplyr::mutate(type=ifelse(type=="exp", "Exponential", "Logistic")) %>%
    dplyr::mutate(best_model=ifelse(best_model=="exp", "Exponential", "Logistic")) %>%
    dplyr::filter(Identity %in% highlight)

  # color_palette = c("firebrick","steelblue"); names(color_palette) = c("Exponential","Logistic")
  color_palette = RColorBrewer::brewer.pal(get_lineages(x) %>% length(), "Dark2")

  if (!show_best) {
    p = rates %>%
      # ggplot(aes(x=Identity, y=rate, ymax=rate, ymin=0, color=ifelse(type==best_model, type, "gainsboro"))) +
      ggplot(aes(x=Identity, group=Lineage, color=Lineage, y=rate, ymax=rate, ymin=0,
                 linetype=type)) +
      geom_linerange(alpha=.8, position=position_dodge2(width=.5)) +
      geom_point(alpha=.8, position=position_dodge2(width=.5)) +

      # facet_wrap(~Lineage) +
      scale_color_manual(values=color_palette) +
      scale_linetype_manual(values=c("solid","dashed") %>% setNames(c("Exponential","Logistic"))) +
      my_ggplot_theme() + theme(panel.grid.major.x=element_blank()) +
      xlab("Clusters") + ylab("Growth rate") + labs(color="Lineage", linetype="")

    return(p)
  }

  dd = rates %>%
    dplyr::filter(type==best_model)

  if (length(dd$best_model %>% unique()) == 1)
    linetps = c("solid") %>% setNames(dd$best_model %>% unique())
  else
    linetps = c("solid","dashed") %>% setNames(c("Logistic","Exponential"))

  p = dd %>%
    ggplot(aes(x=Identity, group=Lineage, color=Lineage, y=rate, ymax=rate, ymin=0, linetype=type)) +
    geom_linerange(alpha=.8, position=position_dodge2(width=.5)) +
    geom_point(alpha=.8, position=position_dodge2(width=.5)) +

    # facet_wrap(~Lineage) +
    scale_color_manual(values=color_palette) +
    scale_linetype_manual(values=linetps) +
    my_ggplot_theme() + theme(panel.grid.major.x=element_blank()) +
    xlab("Clusters") + ylab("Growth rate") + labs(color="Lineage", linetype="")

  return(p)
}
