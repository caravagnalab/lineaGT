# Functions to access elements of a mvnmm object.
get_dataframe = function(x) {
  try(expr = { dataframe = x$dataframe; if (!purrr::is_empty(dataframe)) return(dataframe) }, silent = T)

  py_model = get_model(x)
  return(get_python_dataframe(py_model))
}


get_vaf_dataframe = function(x) {
  return(x$vaf_dataframe)
}


get_viber_clusters = function(x, clusters) {
  if (purrr::is_empty(clusters)) return(get_unique_viber_labels(x))
  vaf = get_vaf_dataframe(x) %>% filter(labels %in% clusters)
  return(vaf$labels_mut %>% unique())
}


get_model = function(x) {
  if ("mvnmm" %in% class(x)) return(x$py_model)
  if ("pylineaGT.mvnmm.MVNMixtureModel" %in% class(x)) return(x)
}

