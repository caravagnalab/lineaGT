.onLoad <- function(libname, pkgname) {
  # unset reticulate python environment, for more details, see:
  # https://github.com/greta-dev/greta/issues/444
  Sys.unsetenv("RETICULATE_PYTHON")

  if (have_lineagt_conda_env()) {
    use_lineagt_conda_env()
  }

}
