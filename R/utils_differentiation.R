get_mrca_df = function(x, highlight, edges, tps=c(), time_spec=F, clonal=F, thr=1, time_dep=T) {
  # time_dep : if TRUE -> assume dependency among timepoints, so a cluster is assigned to a branch based on observations across ALL timepoints
  #            if FALSE -> counts the number of clones in each branch at each timepoint, independently on the previous/next timepoits

  if (length(tps)>0) time_spec = T

  if (purrr::is_empty(tps)) tps = x %>% get_timepoints()

  if (clonal) {
    fracs = x %>%
      get_mean_long() %>%
      dplyr::filter(labels %in% highlight) %>%
      dplyr::filter(timepoints %in% tps) %>%
      dplyr::rename(cluster=labels) %>%
      dplyr::mutate(cluster=as.character(cluster))

    if (time_dep)
      fracs = fracs %>%
        dplyr::group_by(cluster, lineage) %>%
        dplyr::mutate(is.present=any(mean_cov > thr)) %>%
        dplyr::ungroup()
    else
      fracs = fracs %>%
        dplyr::group_by(cluster, lineage, timepoints) %>%
        dplyr::mutate(is.present=any(mean_cov > thr)) %>%
        dplyr::ungroup()

    fracs = fracs %>%
        dplyr::rename(Identity=lineage, Generation=timepoints, Population=mean_cov) %>%
        convert_tp(get_tp_to_int(x)) %>%
        dplyr::inner_join(edges, by="Identity")

    } else if (have_muts_fit(x)) {
      fracs = x %>%
        get_pop_df() %>%
        dplyr::filter(Parent %in% highlight) %>%
        dplyr::filter(Generation %in% get_tp_to_int(x)[tps]) %>%
        # convert_tp(setNames(names(get_tp_to_int(x)), get_tp_to_int(x))) %>%
        dplyr::select(Identity, Population, Lineage, Generation) %>%
        dplyr::rename(cluster=Identity)

      if (time_dep)
        fracs = fracs %>%
          # is.present stores whether the cluster has been observed in at least one tp
          dplyr::group_by(cluster, Lineage) %>%
          dplyr::mutate(is.present=any(Population > thr)) %>%
          dplyr::ungroup()
      else
        fracs = fracs %>%
          # is.present stores whether the cluster has been observed in a specific timepoint
          dplyr::group_by(cluster, Lineage, Generation) %>%
          dplyr::mutate(is.present=any(Population > thr)) %>%
          dplyr::ungroup()

      fracs = fracs %>%
        dplyr::rename(Identity=Lineage) %>%
        dplyr::inner_join(edges, by="Identity")

      } else fracs = data.frame()

  if (nrow(fracs) == 0)
    return(NULL)

  orig = fracs %>%

    # check for each cluster if is present in all descendants of the parent
    dplyr::rowwise() %>%
    dplyr::mutate(is.present.parent=is_present_desc(cluster, Generation, Parent, edges, fracs)) %>%
    dplyr::ungroup() %>%

    # !! cannot filter now otherwise retrieving the mrca won't work
    # dplyr::filter(is.present, Population > thr) %>%

    # orig stores the node of differentiation tree where the cluster has originated
    # it will stay NA if the it's never observed in the lineage
    dplyr::group_by(Generation) %>%
    dplyr::mutate(orig.node=ifelse( (is.present & !is.present.parent), Identity, NA)) %>%
    dplyr::mutate(orig.node=ifelse( (is.na(orig.node) & (is.present & is.present.parent)), Parent, orig.node )) %>%
    dplyr::ungroup() %>%

    dplyr::select(-dplyr::contains("is.present"))

  if (time_dep)
    orig = orig %>%
      dplyr::group_by(Generation, cluster) %>%
      dplyr::mutate(orig.node=ifelse(Generation==Generation & cluster==cluster,
                                get_mrca_list(unique(cluster), edges,
                                              dplyr::filter(., cluster==unique(cluster),Generation==unique(Generation))),
                                orig.node)) %>%
      dplyr::ungroup() %>%

      dplyr::select(-Identity, -Parent) %>%
      dplyr::rename(Identity=orig.node)
  else
    orig = orig %>%
      dplyr::filter(!is.na(orig.node)) %>%
      dplyr::select(-Identity, -Parent) %>%
      dplyr::rename(Identity=orig.node)

  # mrca.list = lapply(unique(orig$Generation), function(gen)
  #   lapply(unique(orig$cluster),
  #          function(cls) {
  #            mrca = get_mrca_list(cls, edges, orig %>% dplyr::filter(Generation==gen))
  #            orig %>% dplyr::mutate(orig=ifelse(Generation==gen & cluster==cls, mrca, orig))
  #          } ) ) %>%
  #     # setNames(nm=unique(orig$cluster)) %>% unlist() ) %>%
  #   data.frame() %>% tibble::rownames_to_column() %>%
  #   setNames(c("cluster",unique(orig$Generation))) %>%
  #   reshape2::melt(id="cluster", value.name="Identity", variable.name="Generation") %>%
  #   dplyr::filter(!is.na(Identity)) %>%
  #   dplyr::inner_join(orig %>% dplyr::select(cluster,Generation,Population) %>% unique())

  # if (all(mrca.list %>% unlist() %>% unique() %>% is.na()))
  if (nrow(orig)==0)
    return(NULL)

  if (time_spec)
    return(
      orig %>%
        dplyr::filter(Population>=thr) %>%
        dplyr::select(-Population) %>% unique() %>%

        dplyr::inner_join(edges, by="Identity") %>%
        dplyr::rename(mrca.from=Parent, mrca.to=Identity) %>%

        dplyr::group_by(mrca.to, Generation) %>%
        dplyr::mutate(n_clones=length(cluster), cluster=paste(cluster, collapse=", ")) %>%
        dplyr::ungroup() %>%

        unique()
    )

  return(
    orig %>%
      dplyr::filter(Population>=thr) %>%
      dplyr::select(-Generation, -Population) %>% unique() %>%
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
  nodes = orig %>% filter(cluster==cls) %>% filter(!is.na(orig.node)) %>% dplyr::pull(orig.node) %>% unique()
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


is_present_desc = function(cls, tp, parent, edges, fracs) {

  desc = get_desc_list(edges)[[parent]]

  # check for each descendant whether cluster "cluster" is observed in lineage "dd" -> Identity contains the "mrca.to"
  for (dd in desc) {
    if ( dd %in% fracs$Identity &&
         !(fracs %>% dplyr::filter(cluster==cls, Generation==tp, Identity==dd) %>% dplyr::pull(is.present)) )
    return(FALSE)
  }
  return(TRUE)
}


get_desc_list = function(edges, clonal=F) {
  desc = list()

  for (pp in unique(edges$Parent)) {  # parents
    if (clonal && pp != "P" && get_parent(edges, pp) != "P")
      next
    for (cc in unique(edges$Identity))  # children
      if ((pp != cc) & (is_desc_of(edges, cc, pp)))
        desc[[pp]] = c(desc[[pp]], cc)  # if cc descends from pp
  }

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

