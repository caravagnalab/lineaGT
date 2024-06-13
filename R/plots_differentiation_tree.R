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
                                     by_timepoint=F,
                                     clonal=T,
                                     wrap=T,
                                     timepoints=c(),
                                     min_abundance=0) {

  highlight = get_highlight(x, highlight=highlight, mutations=clonal)
  if (purrr::is_empty(timepoints)) timepoints = x %>% get_tp_to_int() %>% names()
  if (length(intersect(timepoints, get_timepoints(x)))==0) return(NULL)

  mrca.df = get_mrca_df(x, edges, highlight, tps=timepoints,
                        time_spec=by_timepoint, thr=min_abundance)

  if (!by_timepoint)
    plots = util_plot_diff(mrca.df, edges) else
    plots = lapply(mrca.df$Generation %>% unique(), function(gg)
      util_plot_diff(mrca.df %>% dplyr::filter(Generation==gg), edges=edges) +
        labs(title=paste0("Timepoint ", gg))
      ) %>%
    setNames(mrca.df$Generation %>% unique())


  return(plots)
}


util_plot_diff = function(mrca.df,
                          edges,
                          cex=1,
                          node_palette=colorRampPalette(RColorBrewer::brewer.pal(n=9, "Set1"))) {
  # if (cls != "") {
  #   mrca.df = mrca.list
  #   if (cls %in% names(mrca.df)) mrca.df = mrca.df[[cls]]
  #   cls = paste0("Cluster ", cls, ", ")
  #   }
  # else mrca.df = mrca.list

  if (is.null(mrca.df)) return(ggplot() + theme_minimal())

  mrca.df = mrca.df %>% dplyr::select(-dplyr::contains("branch")) %>%
    dplyr::mutate(cluster=Identity, Identity=mrca.to) %>%
    dplyr::rename(from=mrca.from, to=mrca.to) %>%
    dplyr::full_join(edges %>%
                       dplyr::rename(from=Parent, to=Identity),
                     by=c("from","to")) %>%
    dplyr::mutate(n_clones=replace(n_clones, is.na(n_clones), 0),
                  Identity=to,
                  cluster=replace(cluster, is.na(cluster), ""))

  gg = igraph::graph_from_data_frame(mrca.df)
  igraph::E(gg)$name = mrca.df$to
  igraph::E(gg)$cluster = mrca.df$cluster

  names_vert = igraph::V(gg)$name
  igraph::V(gg)$cluster = dplyr::left_join(data.frame(Identity=names_vert),
                                           dplyr::select(mrca.df, Identity, cluster) %>% unique(),
                                           by="Identity") %>% dplyr::pull(cluster)

  n_nodes = igraph::V(gg) %>% length()
  col_palette = node_palette(n_nodes) %>% setNames(c(edges$Parent, edges$Identity) %>% unique())

  pp = ggraph::ggraph(gg, layout="tree") +
    ggraph::geom_edge_link(aes(width=n_clones, color=name),
                           end_cap=ggraph::circle(5*cex, "mm"),
                           start_cap=ggraph::circle(5*cex, "mm")) +

    ggraph::geom_node_text(aes(label=name), color="black", vjust=0.4) +

    ggrepel::geom_label_repel(aes(label=cluster, x=x, y=y),
                              na.rm=TRUE,
                              nudge_x=.2,
                              nudge_y=.2,
                              size=2.5*cex) +

    ggraph::scale_edge_color_manual(values=col_palette, guide="none") +
    ggraph::scale_edge_width_continuous(name="") + #, breaks=seq(from=.1,to=1.,length.out=15)) +

    guides(color="none") +
    theme_bw() +
    theme(plot.margin=margin(), panel.grid=element_blank(), panel.border=element_blank(),
          axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank(),
          legend.position="bottom")

  # if ("Generation" %in% colnames(mrca.df))
  #   pp = pp + ggraph::facet_edges(Generation ~ ., as.table=T, ncol=4)


  # tb_adj = tidygraph::as_tbl_graph(edges, node_key="Identity")
  #
  # edges.num = tb_adj %>%
  #   tidygraph::activate(edges) %>%
  #   tibble::as_tibble() %>%
  #   tibble::add_column(edges)
  #
  # tb_adj = tb_adj %>%
  #   tidygraph::activate(nodes) %>%
  #   dplyr::rename(Identity=name) %>%
  #   dplyr::left_join(mrca.df, by=c("Identity")) %>%
  #
  #   tidygraph::activate(edges) %>%
  #   dplyr::left_join(edges.num, by=c("from","to")) %>%
  #   dplyr::left_join(mrca.df, by=c("Identity")) %>%
  #
  #   dplyr::mutate(n_clones=ifelse( is.na(n_clones), 0, n_clones )) %>%
  #   dplyr::select(-dplyr::starts_with("mrca"), -Parent)
  #
  # n_nodes = c(edges$Parent, edges$Identity) %>% unique() %>% length()
  # # edge_range = c(min(mrca.df$n_clones), max(mrca.df$n_clones))
  # edge_range = c(1, max(mrca.df$n_clones)+1)
  # col_palette = node_palette(n_nodes) %>% setNames(c(edges$Parent, edges$Identity) %>% unique())
  #
  # layout = ggraph::create_layout(tb_adj, layout="tree")
  # pp = layout %>%
  #   ggraph::ggraph() +
  #   ggraph::geom_edge_link(aes(width=n_clones, colour=Identity),
  #                          end_cap=ggraph::circle(5 * cex, "mm"),
  #                          start_cap=ggraph::circle(5 * cex, "mm")) +
  #
  #   ggraph::geom_node_text(aes(label=Identity), color="black", vjust=0.4) +
  #
  #   ggrepel::geom_label_repel(aes(label=cluster, x=x, y=y),
  #                             na.rm=TRUE,
  #                             nudge_x=.2,
  #                             nudge_y=.2,
  #                             size=2.5*cex) +
  #
  #   ggraph::scale_edge_color_manual(values=col_palette, guide="none") +
  #   ggraph::scale_edge_width_continuous(range=edge_range, breaks=unique(sort(purrr::discard(layout$n_clones, is.na))),
  #                                       name="Cluster Number") +
  #
  #   theme_void(base_size=9*cex) +
  #   theme(legend.position="right", text=element_text(size=9)) +
  #   guides(color="none") +
  #   # labs(title=paste0(cls, "Timepoints ", paste0(timepoints, collapse=","))) +
  #   coord_fixed(0.4) + theme(plot.margin=margin())

  # if ("Generation" %in% colnames(mrca.df))
  #   pp = pp + ggraph::facet_edges(Generation ~ ., as.table=T, ncol=4)
  # else
  #   pp = pp + ggraph::facet_edges( ~ , as.table=T, ncol=4)

  return(pp)
}


differentiation_tree = function(return.numeric=F) {
  edges.diff = data.frame("Parent"=c("WT","P1","P1","P2","P2"), "Identity"=c("P1","P2","Myeloid","B","T"))
  if (return.numeric) { edges.diff$from = c(1,2,2,3,3); edges.diff$to = c(2,3,4,5,6) }
  return(edges.diff)
}
