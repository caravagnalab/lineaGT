#' Visualize the number of subclones on the differentiation tree
#'
#' @param x add
#' @param edges add
#' @param highlight add
#' @param single_tree add
#' @param wrap add
#'
#' @return A \code{ggplot} object showing the identified subclones in the hematopoietic tree.
#'
#' @export plot_differentiation_tree

plot_differentiation_tree = function(x,
                                     edges=differentiation_tree(),
                                     highlight=c(),
                                     single_tree=T,
                                     wrap=T,
                                     timepoints=c()) {

  highlight = get_highlight(x, highlight=highlight)
  if (purrr::is_empty(timepoints)) timepoints = x %>% get_tp_to_int() %>% names()

  if (length(intersect(timepoints, get_timepoints(x)))==0) return(NULL)

  if (length(highlight)==1) cls = highlight else cls = ""

  if (single_tree || length(highlight)==1)
    return(
      util_plot_diff(mrca.list=get_mrca_df(x, highlight=highlight, edges=edges,
                                           tps=timepoints, time_spec=!single_tree),
                     edges=edges,
                     # timepoints=timepoints,
                     cls=cls)
    )

  mrca.list = lapply(highlight, get_mrca_df, x=x, edges=edges, tps=timepoints) %>% setNames(highlight)

  if (mrca.list %>% unlist() %>% unique() %>% is.null()) return(NULL)

  plots = lapply(names(mrca.list), util_plot_diff, mrca.list=mrca.list, edges=edges, timepoints=timepoints) %>%
    purrr::discard(is.null)

  if (wrap) return(patchwork::wrap_plots(plots, guides="collect"))

  return(plots)
}


util_plot_diff = function(cls,
                          mrca.list,
                          edges,
                          # timepoints,
                          cex=1,
                          node_palette=colorRampPalette(RColorBrewer::brewer.pal(n=9, "Set1"))) {
  if (cls != "") {
    mrca.df = mrca.list
    if (cls %in% names(mrca.df)) mrca.df = mrca.df[[cls]]
    cls = paste0("Cluster ", cls, ", ")
    }
  else mrca.df = mrca.list

  if (is.null(mrca.df)) return(ggplot() + theme_minimal())

  mrca.df = mrca.df %>% dplyr::mutate(Identity=mrca.to)

  tb_adj = tidygraph::as_tbl_graph(edges, node_key="Identity")

  edges.num = tb_adj %>%
    tidygraph::activate(edges) %>%
    tibble::as_tibble() %>%
    tibble::add_column(edges)

  tb_adj = tb_adj %>%
    tidygraph::activate(nodes) %>%
    dplyr::rename(Identity=name) %>%
    dplyr::left_join(mrca.df, by=c("Identity")) %>%

    tidygraph::activate(edges) %>%
    dplyr::left_join(edges.num, by=c("from","to")) %>%
    dplyr::left_join(mrca.df, by=c("Identity")) %>%

    dplyr::mutate(n_clones=ifelse( is.na(n_clones), 0, n_clones )) %>%
    dplyr::select(-dplyr::starts_with("mrca"), -Parent)

  n_nodes = c(edges$Parent, edges$Identity) %>% unique() %>% length()
  # edge_range = c(min(mrca.df$n_clones), max(mrca.df$n_clones))
  edge_range = c(1, max(mrca.df$n_clones)+1)
  col_palette = node_palette(n_nodes) %>% setNames(c(edges$Parent, edges$Identity) %>% unique())

  layout = ggraph::create_layout(tb_adj, layout="tree")
  pp = layout %>%
    ggraph::ggraph() +
    ggraph::geom_edge_link(aes(width=n_clones, colour=Identity),
                           end_cap=ggraph::circle(5 * cex, "mm"),
                           start_cap=ggraph::circle(5 * cex, "mm")) +

    ggraph::geom_node_text(aes(label=Identity), color="black", vjust=0.4) +

    ggrepel::geom_label_repel(aes(label=cluster, x=x, y=y),
                              na.rm=TRUE,
                              nudge_x=.2,
                              nudge_y=.2,
                              size=2.5*cex) +

    ggraph::scale_edge_color_manual(values=col_palette, guide="none") +
    ggraph::scale_edge_width_continuous(range=edge_range, breaks=unique(sort(purrr::discard(layout$n_clones, is.na))),
                                        name="Cluster Number") +

    theme_void(base_size=9*cex) +
    theme(legend.position="right", text=element_text(size=9)) +
    guides(color="none") +
    # labs(title=paste0(cls, "Timepoints ", paste0(timepoints, collapse=","))) +
    coord_fixed(0.4) + theme(plot.margin=margin())

  if ("Generation" %in% colnames(df))
    pp = pp + ggraph::facet_edges(Generation ~ vec_gt, as.table=T, ncol=4)
  else
    pp = pp + ggraph::facet_edges( ~ vec_gt, as.table=T, ncol=4)

  return(pp)
}


differentiation_tree = function(return.numeric=F) {
  edges.diff = data.frame("Parent"=c("WT","P1","P1","P2","P2"), "Identity"=c("P1","P2","Myeloid","B","T"))
  if (return.numeric) { edges.diff$from = c(1,2,2,3,3); edges.diff$to = c(2,3,4,5,6) }
  return(edges.diff)
}
