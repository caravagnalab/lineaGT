# Function to retrieve a dataframe with infos about each cluster on the differentation tree.
#

get_mrca_df = function(x, edges, highlight=c(), tps=c(), time_spec=F, thr=0) {
  # time_dep : if TRUE -> assume dependency among timepoints, so a cluster is assigned to a branch based on observations across ALL timepoints
  #            if FALSE -> counts the number of clones in each branch at each timepoint, independently on the previous/next timepoints

  if (purrr::is_empty(highlight)) highlight = get_highlight(x, mutations=T)
  if (purrr::is_empty(tps)) tps = x %>% get_tp_to_int()
  if (is.character(tps)) tps = get_tp_to_int(x)[tps]  # take the integer values

  # get population sizes and phylogenies (clone tree) for each generation
  fracs = x %>%
    get_pop_df() %>%
    dplyr::select(Identity, Generation, Lineage, Population, Parent) %>%
    dplyr::filter(Identity %in% highlight) %>%
    dplyr::filter(Generation %in% tps) %>%
    dplyr::mutate(Identity=as.character(Identity)) %>%
    dplyr::arrange(Identity, Lineage, Generation) %>%

    tibble::add_column(mrca.to=NA, next.tp.mrca=NA) %>%
    dplyr::filter(Population>thr)

  if (time_spec)
    # retrieves info about location of each clone on the differentiation tree
    fracs = get_time_spec_mrca(fracs, tps, edges)
  else
    fracs = fracs %>%
      dplyr::group_by(Identity) %>%
      dplyr::mutate(mrca.to=lowest_common_acestor(Lineage, edges)) %>%
      dplyr::ungroup()

  fracs = fracs %>%
    dplyr::inner_join(edges %>% dplyr::rename(mrca.to=Identity, mrca.from=Parent), by="mrca.to") %>%
    dplyr::mutate(branch=paste0(mrca.from,"->",mrca.to)) %>%
    dplyr::mutate(Identity=factor(Identity, levels=highlight)) %>%
    dplyr::arrange(Identity)


  if (time_spec)
    return(
      fracs %>%
        dplyr::select(Identity, Generation, dplyr::contains(c("branch","mrca"))) %>% unique() %>%
        dplyr::group_by(Generation, branch) %>%
        dplyr::reframe(n_clones=sum(!grepl(pattern="S", x=unique(Identity))),
                       n_subclones=sum(grepl(pattern="S", x=unique(Identity))),
                       Identity=toString(unique(Identity)),
                       mrca.from=mrca.from, mrca.to=mrca.to) %>% unique() %>%
        dplyr::select(dplyr::contains("from"), dplyr::contains("to"), dplyr::everything())
    )

  return(
    fracs %>%
      dplyr::select(Identity, dplyr::contains(c("branch","mrca"))) %>% unique() %>%
      dplyr::group_by(branch) %>%
      dplyr::reframe(n_clones=sum(!grepl(pattern="S", x=unique(Identity))),
                     n_subclones=sum(grepl(pattern="S", x=unique(Identity))),
                     Identity=toString(unique(Identity)),
                     mrca.from=mrca.from, mrca.to=mrca.to) %>% unique() %>%
      dplyr::select(dplyr::contains("from"), dplyr::contains("to"), dplyr::everything())
  )
}


get_time_spec_mrca = function(fracs, tps, edges) {
  # start from the last generation
  sorted_gens = sort(as.array(tps), decreasing=T)

  for (gg in 1:length(sorted_gens)) {
    fracs = fracs %>%

      # compute the "location" of each cluster in each timepoint based on population size
      dplyr::group_by(Identity, Generation) %>%
      dplyr::mutate(mrca.to=ifelse(Generation==sorted_gens[gg],
                                   lowest_common_acestor(c(Lineage, unique(next.tp.mrca)), edges),
                                   mrca.to)) %>%
      dplyr::ungroup() %>%

      # annotate the next analysed timepoint (i.e., previous time) with the current mrca
      get_mrca_next_tp(sorted_gens, gg)
  }

  return(fracs %>% dplyr::select(-next.tp.mrca))
}


# get_mrca_list = function(cls, edges, orig) {
#   # cls is a subclone
#   # edges is the edges dataframe
#   # orig is a dataframe reporting where each cluster has been observed
#
#   # the function's purpose is to return the MRCA of "cls"
#   nodes = orig %>% filter(cluster==cls) %>% filter(!is.na(orig.node)) %>% dplyr::pull(orig.node) %>% unique()
#   root = get_root(edges)
#
#   if (length(unique(nodes)) == 1) return(nodes %>% unique())
#
#   mrca = nodes[1]
#
#   for (n1 in nodes)
#     for (n2 in nodes) {
#       if (n1 == n2) next
#       tmp = get_mrca(edges, n1, n2)
#       mrca = get_mrca(edges, mrca, tmp)
#       if (mrca == root) return(root)
#     }
#   return(mrca)
# }


# get_mrca = function(edges, n1, n2) {
#   root = get_root(edges)
#   if (n1==n2 && n1==root) return(root)
#   if (n1==n2) return(n1)
#
#   p1 = get_parent(edges, n1)
#   p2 = get_parent(edges, n2)
#
#   if (is_desc_of(edges,n1,n2)) return(n2)  # n1 descends from n2
#   if (is_desc_of(edges,n2,n1)) return(n1)  # n2 descends from n1
#   if (is_desc_of(edges,n1,p2)) return(p2)  # n1 descends from p2
#   if (is_desc_of(edges,n2,p1)) return(p1)  # n2 descends from p1
#
#   return(get_mrca(edges, p1, p2))  # check for their parents
# }


# is_present_desc = function(cls, tp, parent, edges, fracs) {
#
#   desc = get_desc_list(edges)[[parent]]
#
#   # check for each descendant whether cluster "cluster" is observed in lineage "dd" -> Identity contains the "mrca.to"
#   for (dd in desc) {
#     if ( dd %in% fracs$Identity &&
#          !(fracs %>% dplyr::filter(cluster==cls, Generation==tp, Identity==dd) %>% dplyr::pull(is.present)) )
#     return(FALSE)
#   }
#   return(TRUE)
# }


get_desc_list = function(edges, clonal=F, leaves=F) {
  desc = list()

  for (pp in unique(edges$Parent)) {  # parents
    if (clonal && pp != "P" && get_parent(edges, pp) != "P")
      next
    for (cc in unique(edges$Identity))  # children
      if ((pp != cc) & (is_desc_of(edges, cc, pp)))
        desc[[pp]] = c(desc[[pp]], cc)  # if cc descends from pp
  }

  if (leaves) for (ll in setdiff(edges$Identity, edges$Parent)) desc[[ll]] = ll

  return(desc)
}


# compute_n_clones = function(edges, mrca.df, id) {
#   if (id %in% mrca.df$mrca.to) n_clones = mrca.df[mrca.df$mrca.to==id,] %>% dplyr::pull(n_clones)
#   else n_clones = 0
#
#   for (nn in mrca.df$mrca.to)
#     if (is_desc_of(edges, id, nn)) n_clones = n_clones + mrca.df[mrca.df$mrca.to==nn,] %>% dplyr::pull(n_clones)
#
#   return(n_clones)
# }


edges_to_matrix = function(edges) {
  vars = unique(unlist(edges))
  matrix = matrix(0, nrow = length(vars), ncol = length(vars))
  colnames(matrix) = vars
  rownames(matrix) = vars

  if(nrow(edges) == 0) return(matrix)

  for(j in 1:nrow(edges)) matrix[edges[j, "Parent"], edges[j, "Identity"]] = 1

  return(matrix)
}


# matrix_to_edges = function(matr) {
#   dfedges = data.frame(stringsAsFactors = F)
#   for(i in 1:nrow(matr))
#     for(j in 1:ncol(matr))
#       if(matr[i,j] == 1)
#         dfedges = rbind(dfedges, data.frame(Parent=rownames(matr)[i], Identity=colnames(matr)[j], stringsAsFactors=FALSE))
#
#   return(dfedges)
# }


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


lowest_common_acestor = function(nodes, edges) {
  # nodes is a list of nodes for which we want to find the mrca
  nodes = nodes[!is.na(nodes)]

  # get for each node in "edges", the list of descendants
  desc_list = get_desc_list(edges, leaves=T)

  return(
    lapply(names(desc_list), function(x) {
      # if all nodes are descendants of "x", then returns the number of descendants, NULL otherwise
      if (all(nodes %in% c(x, desc_list[[x]]))) return(length(desc_list[[x]]))
      }) %>%
      setNames(names(desc_list)) %>% purrr::discard(is.null) %>%
      # keep the node with less descendants -> mrca
      which.min() %>% names
  )
}



get_mrca_next_tp = function(fracs, tps, tp.idx) {
  # is.last.tp tells me if "tp" is the last one, if so, the new column will be NA
  # tps = list of timepoints

  # if we are dealing with the first timepoint
  if (tp.idx == length(tps))
    return(fracs)

  # save the state at the current timepoint
  fracs.curr = fracs %>%
    dplyr::select(Identity, Generation, mrca.to) %>% unique() %>%
    dplyr::filter(Generation==tps[tp.idx]) %>%

    dplyr::rowwise() %>%
    dplyr::mutate(Generation=get_next_tp(fracs, tps[tp.idx], Identity)) %>%
    dplyr::ungroup() %>%

    dplyr::rename(next.tp.mrca2=mrca.to)

  return(
    fracs %>%
      dplyr::ungroup() %>%
      dplyr::left_join(fracs.curr, by=c("Identity", "Generation")) %>%
      dplyr::mutate(next.tp.mrca=ifelse(is.na(next.tp.mrca), next.tp.mrca2, next.tp.mrca)) %>%
      dplyr::select(-next.tp.mrca2)
  )
}

get_next_tp = function(fracs, gen, idd) {
  # print(c(idd, gen))
  # print(fracs %>%
  #         dplyr::filter(Identity==idd) %>%
  #         dplyr::filter(Generation < gen))

  next.gens = fracs %>%
    dplyr::filter(Identity==idd) %>%
    dplyr::filter(Generation < gen) %>%
    dplyr::pull(Generation)

  if (length(next.gens) == 0)
    return(NA)

  return(next.gens %>% max())
}

fix_missing_tps = function(fracs, tps, tp.idx) {
  tmp = fracs %>%
    dplyr::rowwise() %>%
    dplyr::mutate(is.missing=is.na(next.tp.mrca) && Generation==tps[tp.idx+1]) %>%
    dplyr::group_by(Identity) %>%
    dplyr::filter(any(is.missing)) %>%
    dplyr::ungroup()

  tmp %>%
    dplyr::group_by(Identity) %>%
    dplyr::mutate(next.tp=ifelse(Generation==tps[tp.idx+1],
                                 get_next_tp(Generation),
                                 ))

}






