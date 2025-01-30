#' Muller plot
#'
#' @description Function to visualize the mullerplot for the fitted object.
#'
#' @param x a mvnmm object.
#' @param which string among \code{"frac","pop","fitness"} determining whether to plot the coverage
#' normalized in \code{[0,1]}, as absolute clone abundance, or with each clone colored by the growth rate,
#' computed assuming a exponential growth.
#' @param highlight a vector of clusters IDs to highlight in the plot.
#' @param min_frac min_frac numeric value in \code{[0,1]} representing the minimum abundance to highlight a clone.
#' @param timepoints_to_int a list to map each \code{timepoint} value to an integer.
#' @param mutations Boolean. If set to \code{TRUE}, also the clusters of mutations will be visualized.
#' @param single_clone Boolean. If \code{mutations} and \code{single_clone} are set to \code{TRUE}, only the clones
#' reported in \code{highlight} and the respective subclones will be visualised.
#' @param rm_mixt remove clusters estimated as polyclonal
#' @param tree_score add
#' @param legend.pos add
#'
#' @examples
#' if (FALSE) plot_mullerplot(x, wrap=T)
#'
#' @import ggplot2
#' @import ggmuller
#' @importFrom patchwork wrap_plots
#'
#' @export plot_mullerplot

plot_mullerplot = function(x,
                           which="frac",
                           highlight=c(),
                           min_frac=0,
                           min_abundance=0,
                           estimate_npops=FALSE,
                           vcn=NULL,
                           rm_mixt=FALSE,
                           timepoints_to_int=c(),
                           mutations=F,
                           single_clone=T,
                           tree_score=1,
                           legend.pos="right") {

  timepoints_to_int = map_timepoints_int(x, timepoints_to_int)

  highlight.cov = get_highlight(x, min_frac=min_frac, highlight=highlight)

  if (rm_mixt) {
    n_pops = estimate_n_pops(x, highlight=highlight.cov, vcn=vcn)

    if (length(names(n_pops[n_pops==1])) == length(highlight.cov))
      cli::cli_alert_info("All the estimated clusters represent a single population. No clusters are removed.")

    if (length(names(n_pops[n_pops==1])) < length(highlight.cov))
      cli::cli_alert_info("Clusters {.field {names(n_pops[n_pops>1])}} represent more than one population. They will be removed from the visualisation.")

    highlight.cov = intersect(highlight.cov, names(n_pops[n_pops==1]))
  }

  highlight = get_highlight(x, min_frac=min_frac, min_abundance=min_abundance, highlight=highlight.cov, mutations=mutations)
  color_palette = highlight_palette(x, highlight)
  lvls = c("P", get_unique_muts_labels(x), get_unique_labels(x))

  pop_df = get_muller_pop(x,
                          mutations=mutations,
                          timepoints_to_int=timepoints_to_int,
                          highlight=highlight.cov,
                          single_clone=single_clone,
                          estimate_npops=estimate_npops,
                          vcn=vcn) %>%
    dplyr::select(-Population, -Frequency, -Parent, -theta_binom, -dplyr::contains("Pop.subcl")) %>%
    dplyr::rename(Population=Pop.plot) %>%
    dplyr::arrange(Identity, Generation, Lineage)



  if (estimate_npops)
    pop_df = pop_df %>%
      dplyr::rename(Population.orig=Population) %>%
      dplyr::rename(Population=Population.corr)

  edges_df = get_muller_edges(x,
                              mutations=mutations,
                              tree_score=tree_score) %>%
    dplyr::arrange(Parent)


  if (single_clone && mutations) {
    pop_df = pop_df %>% filter_muller_df(highlight=highlight.cov)
    edges_df = edges_df %>% filter_muller_df(highlight=highlight.cov)
  }

  timepoints = x %>% get_dimensions()
  lineages = x %>% get_lineages()

  mullerdf = data.frame()
  for (ll in lineages) {
    tp = timepoints[grep(pattern=ll, x=timepoints)]
    if (length(tp) != 0) {
      pop_ll = pop_df %>% dplyr::filter(Lineage==ll)
      mullerdf_ll = ggmuller::get_Muller_df(edges_df, pop_ll) %>%
        dplyr::mutate(Lineage=ll)

      mullerdf = rbind(mullerdf, mullerdf_ll)
    }
  }

  return(
    mullerplot_util(x, mullerdf %>% filter(Generation %in% (timepoints_to_int %>% unlist())),
                    which=which,
                    highlight=stringr::str_sort(highlight, numeric=TRUE),
                    color_palette=color_palette,
                    legend.pos=legend.pos,
                    estimate_npops=estimate_npops,
                    vcn=vcn)
  )
}


mullerplot_util = function(x, mullerdf, which, color_palette, highlight, legend.pos="right", estimate_npops=F, vcn=NULL) {
  if (!which %in% c("fitness", "frac", "pop")) {
    cli::format_warning('`which` must be one among "frac","pop","fitness". Using `which`="frac".')
    which = "frac"
  }

  if (which == "fitness") {
    mullerdf = x %>% get_growth_rates() %>%
      dplyr::filter(type==best_model) %>%
      dplyr::select(Lineage, Identity, rate) %>%
      dplyr::full_join(mullerdf, by=c("Lineage","Identity")) %>%
      dplyr::mutate(rate=replace(rate, is.na(rate), 0))
    exp_limits = c(min(mullerdf$rate), max(mullerdf$rate))

    return(
      mullerdf %>%
        ggplot() +
        geom_area(aes_string(x="Generation", y="Frequency", group="Group_id", fill="rate")) +
        geom_vline(xintercept=mullerdf$Generation %>% unique(), linetype="dashed") +
        guides(linetype="none", color="none") +
        facet_wrap(~Lineage, nrow=1) +
        xlab("Time") +
        labs(fill="Exp rate") +
        my_ggplot_theme(legend.pos=legend.pos) +
        scale_fill_gradient2(mid="white", low="blue", high="red",
                             limits=exp_limits, na.value="#FFFFFF00")
    )
  }

  if (which == "frac")
    y = "Frequency"
  else if (which == "pop") {
    y = "Population"
    mullerdf = mullerdf %>% pop_df_add_empty()
  }


  fillname = "Clusters"
  if (estimate_npops) {
    pp = x %>% estimate_n_pops(vcn=vcn)
    names(color_palette) = sapply(names(color_palette),
                                  function(a) {
                                    if (a %in% names(pp)) return(paste0(a, " - ", pp[[a]]))
                                    return(a)
                                  } ) %>%
      setNames(NULL)

    highlight = sapply(highlight, function(a) {
      if (a %in% names(pp)) return(paste0(a, " - ", pp[[a]]))
      return(a)
      } ) %>% setNames(NULL)

    fillname = "Clusters - # pops"

    if (!"n_pops" %in% colnames(mullerdf))
      mullerdf = data.frame(pp) %>%
        tibble::rownames_to_column(var="Identity") %>%
        dplyr::inner_join(mullerdf, by="Identity")

    mullerdf = mullerdf %>%
      dplyr::rowwise() %>%
      dplyr::mutate(n_pops=as.character(n_pops),
                    Identity=replace(Identity, (!is.na(n_pops) && n_pops!="0"), paste0(Identity, " - ", n_pops))) %>%
      dplyr::ungroup()
  }

  id_levels = mullerdf %>%
    dplyr::pull(Group_id) %>% levels() %>%
    purrr::discard(function(i) i %in% c("___special_empty","___special_emptya","P","Pa")) %>%
    purrr::discard(function(i) grepl(".S", i)) %>%
    stringr::str_remove_all("a") %>% unique()

  mullerdf_lines = mullerdf %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Identity=dplyr::case_when(
      grepl(".S", Identity) ~ strsplit(Identity, split="[.]")[[1]][1],
      .default=Identity
    )) %>%
    dplyr::mutate(Identity=factor(Identity, levels=id_levels)) %>%
    dplyr::group_by(Generation, Identity, Lineage) %>%
    dplyr::summarise(Population=sum(Population),
                     Frequency=sum(Frequency)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate()

  return(
    mullerdf %>%
      ggplot() +
      geom_area(aes_string(x="Generation",
                           y=y,
                           group="Group_id",
                           fill="Identity"), alpha=1, lwd=0) +
      geom_line(aes_string(x="Generation", y=y, group="Identity"), position="stack",
                inherit.aes=FALSE, data=mullerdf_lines, color="black", linewidth=0.4) +
      geom_vline(xintercept=mullerdf$Generation %>% unique(), linetype="dashed") +
      guides(linetype="none", color="none") +
      facet_wrap(~Lineage, nrow=1) +
      scale_fill_manual(name=fillname, values=color_palette, na.value="white", breaks=highlight) +
      scale_color_manual(name=fillname, values=color_palette, na.value="white", breaks=highlight) +
      xlab("Time") +
      my_ggplot_theme(legend.pos=legend.pos) +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  )
}

