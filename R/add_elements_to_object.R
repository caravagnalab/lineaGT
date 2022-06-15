add_viber_run = function(x, viber_run, label="") {
  if (label == "") x$viber_run = viber_run
  else x[[paste("viber_run", label, sep=".")]] = viber_run
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

