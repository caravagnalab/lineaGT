.onLoad <- function(libname, pkgname) {
  # reticulate::configure_environment(pkgname, force=T)
  reticulate::install_miniconda(path=reticulate::miniconda_path(), update=TRUE, force=TRUE)
  reticulate::py_install("pylineaGT", pip=T)
  print(reticulate::py_available())
}
