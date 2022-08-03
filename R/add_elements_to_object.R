## Set of functions used to add to the mvnmm object

add_muts_fit = function(x, x.muts, label="") {
  if (label == "") x$x.muts = x.muts
  else x[[paste("x.muts", label, sep=".")]] = x.muts
  return(x)
}


add_color_palette = function(x, color_palette, label="") {
  if (label == "") x$color_palette = color_palette
  else x[[paste("color_palette", label, sep=".")]] = color_palette
  return(x)
}


add_phylo = function(x, trees, label="") {
  if (label == "") x$trees = trees
  else x[[paste("trees", label, sep=".")]] = trees
  return(x)
}


add_vaf = function(x, vaf.df, label="") {
  if (label == "") x$vaf.dataframe = vaf.df
  else x[[paste("vaf.dataframe", label, sep=".")]] = vaf.df
  return(x)
}


# add_tp_int = function(x, tp_int) {
#
# }
