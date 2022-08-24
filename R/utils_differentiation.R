get_mrca_df = function(x, highlight, edges, tps=c()) {

  if (purrr::is_empty(tps)) tps = x %>% get_timepoints()

  # edges.diff = differentiation_tree(return.numeric=T)
  fracs = x %>%
    get_vaf_dataframe() %>%
    dplyr::filter(labels %in% highlight) %>%
    dplyr::filter(timepoints %in% tps) %>%
    dplyr::select(labels_mut, theta_binom, lineage, timepoints) %>%
    dplyr::rename(cluster=labels_mut) %>%

    # is.present stores whether the cluster has been observed in the cluster
    dplyr::group_by(cluster, lineage) %>%
    dplyr::summarise(is.present=any(theta_binom > 0), .groups="keep") %>%
    dplyr::ungroup() %>%
    dplyr::rename(Identity=lineage) %>%
    dplyr::inner_join(edges, by="Identity")

  if (nrow(fracs) == 0)
    return(NULL)

  orig = fracs %>%
    dplyr::rowwise() %>%
    dplyr::mutate(is.present.parent=is_present_desc(cluster, Parent, edges, fracs) ) %>%  # if is present in all descendants of the parent
    dplyr::ungroup() %>%

    # orig stores the node of differentiation tree where the cluster has originated
    # it will stay NA if the it's never observed in the lineage
    dplyr::mutate(orig=ifelse( (is.present & !is.present.parent), Identity, NA)) %>%
    dplyr::mutate(orig=ifelse( (is.na(orig) & (is.present & is.present.parent)), Parent, orig ))

  mrca.list = lapply(unique(orig$cluster), get_mrca_list, edges=edges, orig=orig) %>%
    setNames(nm=unique(orig$cluster))

  if (all(mrca.list %>% unlist() %>% unique() %>% is.na()))
    return(NULL)

  return(
    mrca.list %>%
      data.frame() %>% t() %>% as.data.frame() %>%
      tibble::rownames_to_column() %>% setNames(c("cluster","Identity")) %>%
      dplyr::inner_join(edges, by="Identity") %>%
      dplyr::rename(mrca.from=Parent, mrca.to=Identity) %>%

      dplyr::group_by(mrca.to) %>%
      dplyr::mutate(n_clones=length(cluster), cluster=paste(cluster, collapse=", ")) %>%
      dplyr::ungroup() %>%

      unique()
  )
}


get_mrca_list = function(cls, edges, orig) {
  # cls is a subclone
  # edges is the edges dataframe
  # orig is a dataframe reporting where each cluster has been observed

  # the function's purpose is to return the MRCA of "cls"
  nodes = orig %>% filter(cluster==cls) %>% filter(!is.na(orig)) %>% dplyr::pull(orig)
  root = get_root(edges)

  if (length(unique(nodes)) == 1) return(nodes %>% unique())

  mrca = nodes[1]

  for (n1 in nodes)
    for (n2 in nodes) {
      if (n1 == n2) next
      tmp = get_mrca(edges, n1, n2)
      mrca = get_mrca(edges, mrca, tmp)
      if (mrca == root) return(root)
    }
  return(mrca)
}


get_mrca = function(edges, n1, n2) {
  root = get_root(edges)
  if (n1==n2 && n1==root) return(root)
  if (n1==n2) return(n1)

  p1 = get_parent(edges, n1)
  p2 = get_parent(edges, n2)

  if (is_desc_of(edges,n1,n2)) return(n2)  # n1 descends from n2
  if (is_desc_of(edges,n2,n1)) return(n1)  # n2 descends from n1
  if (is_desc_of(edges,n1,p2)) return(p2)  # n1 descends from p2
  if (is_desc_of(edges,n2,p1)) return(p1)  # n2 descends from p1

  return(get_mrca(edges, p1, p2))  # check for their parents
}


is_present_desc = function(node, parent, edges, fracs) {
  desc = get_desc_list(edges)[[parent]]

  # check for each descendant whether cluster "node" is observed in lineage "dd"
  for (dd in desc) {
    if ( dd %in% fracs$Identity &&
         !(fracs %>% dplyr::filter(cluster==node, Identity==dd) %>% dplyr::pull(is.present)) )
      return(FALSE)
  }
  return(TRUE)
}


get_desc_list = function(edges) {
  desc = list()

  for (pp in unique(edges$Parent))  # parents
    for (cc in unique(edges$Identity))  # children
      if ((pp != cc) & (is_desc_of(edges, cc, pp)))
        desc[[pp]] = c(desc[[pp]], cc)  # if cc descends from pp

  return(desc)
}


compute_n_clones = function(edges, mrca.df, id) {
  if (id %in% mrca.df$mrca.to) n_clones = mrca.df[mrca.df$mrca.to==id,] %>% dplyr::pull(n_clones)
  else n_clones = 0

  for (nn in mrca.df$mrca.to)
    if (is_desc_of(edges, id, nn)) n_clones = n_clones + mrca.df[mrca.df$mrca.to==nn,] %>% dplyr::pull(n_clones)

  return(n_clones)
}


edges_to_matrix = function(edges) {
  vars = unique(unlist(edges))
  matrix = matrix(0, nrow = length(vars), ncol = length(vars))
  colnames(matrix) = vars
  rownames(matrix) = vars

  if(nrow(edges) == 0) return(matrix)

  for(j in 1:nrow(edges)) matrix[edges[j, "Parent"], edges[j, "Identity"]] = 1

  return(matrix)
}


matrix_to_edges = function(matr) {
  dfedges = data.frame(stringsAsFactors = F)
  for(i in 1:nrow(matr))
    for(j in 1:ncol(matr))
      if(matr[i,j] == 1)
        dfedges = rbind(dfedges, data.frame(Parent=rownames(matr)[i], Identity=colnames(matr)[j], stringsAsFactors=FALSE))

  return(dfedges)
}


get_parent = function(edges, node) {
  return(
    edges %>%
      dplyr::filter(Identity==node) %>%
      dplyr::pull(Parent)
  )
}


is_desc_of = function(edges, desc, anc) {
  root = get_root(edges)

  if (anc == root) return(TRUE)
  if (desc == root) return(FALSE)

  par = get_parent(edges, desc)
  if (par == anc) return(TRUE)

  return(is_desc_of(edges, par, anc))
}


get_root = function(edges) {
  return(setdiff(edges$Parent, edges$Identity))
}

