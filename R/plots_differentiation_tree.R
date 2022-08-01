#' Visualize the number of subclones on the differentiatio tree
#'
#' @param x add
#' @param label add
#' @param label.on.node add
#' @param cex add
#'
#' @return
#' @export plot_differentiation_tree

plot_differentiation_tree = function(x, label="", label.on.node=T, cex=1) {
  edges.diff = differentiation_tree()
  mrca.df = get_mrca_df(x, label=label)
  desc = get_desc_list(edges.diff)

  if (label.on.node) mrca.df = mrca.df %>% dplyr::mutate(Identity=mrca.to)
  else mrca.df = mrca.df %>% dplyr::mutate(Identity=mrca.from)

  tb_adj = tidygraph::as_tbl_graph(edges.diff) %>%
    tidygraph::activate(nodes) %>%
    dplyr::rename(Identity=name) %>%
    dplyr::left_join(mrca.df, by=c("Identity")) %>%
    dplyr::group_by(Identity) %>%
    dplyr::mutate(n_clones=ifelse( is.na(n_clones), compute_n_clones(edges.diff,mrca.df,id=Identity), n_clones )) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster=ifelse( cluster=="NA", "", cluster )) %>%
    dplyr::select(-from, -to, -dplyr::contains("mrca")) %>%

    tidygraph::activate(edges) %>%
    dplyr::left_join(mrca.df, by=c("from","to")) %>%
    dplyr::mutate(cluster=ifelse( is.na(cluster), "", cluster )) %>%
    dplyr::select(-Identity, -dplyr::starts_with("mrca"), -n_clones)

  # n_nodes = c(edges.diff$Parent, edges.diff$Identity) %>% unique() %>% length()
  # nodes_names = c(edges.diff$Parent, edges.diff$Identity) %>% unique()
  # node_palette = RColorBrewer::brewer.pal(n=n_nodes, "Set1") %>% setNames(nodes_names)

  layout = ggraph::create_layout(tb_adj, layout = "tree")
  if (label.on.node)
    pp = layout %>%
      ggraph::ggraph() +
      ggraph::geom_edge_link(arrow=arrow(length=unit(2 * cex, "mm")),
                             end_cap=ggraph::circle(5 * cex, "mm"),
                             start_cap=ggraph::circle(5 * cex, "mm")) +
      # if name on the nodes, we need to use the mrca.to nodes
      ggrepel::geom_label_repel(aes(label=cluster, x=x, y=y), na.rm=TRUE, nudge_x=.3, nudge_y=.3, size=2.5*cex)

  else
    pp = layout %>%
      ggraph::ggraph() +
      ggraph::geom_edge_link(aes(label=cluster),
                             arrow=arrow(length=unit(2 * cex, "mm")),
                             end_cap=ggraph::circle(5 * cex, "mm"),
                             start_cap=ggraph::circle(5 * cex, "mm"),
                             angle_calc="none",
                             label_dodge = unit(2.5, 'mm'))

  pp = pp +
    # ggrepel::geom_label_repel(aes(label=n_clones, x=x, y=y), na.rm=TRUE, nudge_x=.3, nudge_y=.3, size=2.5*cex) +
    ggraph::geom_node_point(aes(size=n_clones, color=Identity), alpha=.7) +
    ggraph::geom_node_text(aes(label=Identity), colour="black", vjust=0.4) +
    coord_cartesian(clip="off") +
    # scale_color_manual(values=node_palette) +
    scale_size(range=c(3, 10) * cex) +
    guides(color=FALSE, size=guide_legend("Clone size", nrow=1)) +
    theme_void(base_size=10*cex) +
    theme(legend.position="bottom", text=element_text(size=9))

  return(pp)
}


differentiation_tree = function(return.numeric=F) {
  edges.diff = data.frame("Parent"=c("GL","P1","P1","P2","P2"), "Identity"=c("P1","P2","Myeloid","B","T"))
  if (return.numeric) { edges.diff$from = c(1,2,2,3,3); edges.diff$to = c(2,3,4,5,6) }
  return(edges.diff)
}


get_mrca = function(nodes, edges) {
  if (length(nodes) == 1) return(nodes)
  if (length(unique(nodes)) == 1) return(nodes %>% unique())

  mrca = nodes[1]; found = F

  for (n1 in nodes)
    for (n2 in nodes)
      if ((n1 != n2) & (is_desc_of(edges, n1, n2))) {
        mrca = n2; found = T  # is the ancestor
      }

  if (!found) return(setdiff(edges$Parent, edges$Identity))
  return(mrca)
}
