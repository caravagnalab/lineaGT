## Set of functions used to add to the mvnmm object

add_muts_fit = function(x, x.muts, label="") {
  x$x.muts = x.muts
  return(x)
}


add_color_palette = function(x, color.palette, label="") {
  x$color.palette = color.palette
  return(x)
}


add_phylo = function(x, x.trees, label="") {
  x$x.trees = x.trees
  return(x)
}


add_vaf = function(x, vaf.df, label="") {
  x$vaf.dataframe = vaf.df
  return(x)
}


add_growth_rates = function(x, rates.df, label="") {
  x$growth.rates = rates.df
  return(x)
}


add_tp_int = function(x, tp.to.int) {
 if ("tp.to.int" %in% names(x)) return(x)
 x$tp.to.int = tp.to.int
 return(x)
}
