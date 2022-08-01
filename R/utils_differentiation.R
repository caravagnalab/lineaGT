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
  root = setdiff(edges$Parent, edges$Identity)

  if (anc == root) return(TRUE)
  if (desc == root) return(FALSE)

  par = get_parent(edges, desc)
  if (par == anc) return(TRUE)

  return(is_desc_of(edges, par, anc))
}


get_mrca_df_single_clone = function(mut, cloneID="") {
  edges.diff = differentiation_tree(return.numeric=T)

  ccf = get_binomial_theta_cluster(mut, cloneID) %>%
    dplyr::rename(cluster=labels_mut) %>%
    dplyr::group_by(cluster, lineage) %>%
    dplyr::summarise(is.present=any(theta>0), .groups="keep") %>%
    dplyr::rename(Identity=lineage) %>%
    dplyr::inner_join(edges.diff, by="Identity") %>%
    dplyr::ungroup()

  orig = ccf %>%
    dplyr::group_by(Parent, cluster) %>%
    dplyr::mutate(is.present.parent=all(is.present)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(orig=ifelse( (is.present & !is.present.parent), Identity, NA)) %>%
    dplyr::mutate(orig=ifelse( (is.na(orig) & (is.present & is.present.parent)), Parent, orig ))

  mrca.list = list()
  for (cl in orig$cluster %>% unique()) {
    nodes = orig %>% filter(cluster==cl) %>% filter(!is.na(orig)) %>% dplyr::pull(orig)
    mrca.list[[cl]] = get_mrca(nodes, edges.diff)
  }

  return(
    mrca.list %>%
      data.frame() %>% t() %>% as.data.frame() %>%
      rownames_to_column() %>% setNames(c("cluster","Identity")) %>%
      dplyr::inner_join(edges.diff, by="Identity")
  )
}


get_mrca_df = function(x, label="") {
  muts = x %>% get_muts_fit()

  mrca.df = data.frame()
  for (cl in names(muts)) {
    mut = muts[[cl]]
    if (purrr::is_empty(mut)) next

    if (purrr::is_empty(mrca.df)) mrca.df = get_mrca_df_single_clone(mut, cloneID=cl)
    else mrca.df = mrca.df %>% dplyr::add_row( get_mrca_df_single_clone(mut, cloneID=cl) )
  }

  return(
    mrca.df %>%
      dplyr::rename(mrca.from=Parent, mrca.to=Identity) %>%
      group_by(mrca.to) %>%
      dplyr::mutate(n_clones=length(cluster), cluster=paste(cluster, collapse=", ")) %>%
      ungroup() %>%
      unique()
  )
}


get_desc_list = function(edges) {
  desc = list()

  for (pp in unique(edges$Parent))  # parents
    for (cc in unique(edges$Identity))  # children
      if ((pp != cc) & (is_desc_of(edges, cc, pp)))  desc[[pp]] = c(desc[[pp]], cc)  # if cc descends from pp

  return(desc)
}


compute_n_clones = function(edges, mrca.df, id) {
  if (id %in% mrca.df$mrca.to) n_clones = mrca.df[mrca.df$mrca.to==id,] %>% dplyr::pull(n_clones)
  else n_clones = 0

  for (nn in mrca.df$mrca.to)
    if (is_desc_of(edges, id, nn)) n_clones = n_clones + mrca.df[mrca.df$mrca.to==nn,] %>% dplyr::pull(n_clones)

  return(n_clones)
}

