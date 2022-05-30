.onLoad <- function(libname, pkgname) {
  # unset reticulate python environment, for more details, see:
  # https://github.com/greta-dev/greta/issues/444
  Sys.unsetenv("RETICULATE_PYTHON")

  # if (have_lineagt_conda_env()) {
  #   cat("The conda environment 'lineagt_env' is present and will be loaded!")
  #   use_conda_env()
  # } else {
  configure_environment()
  # }

}
