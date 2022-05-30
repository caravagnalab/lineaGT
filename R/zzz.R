.onLoad <- function(libname, pkgname) {
  # unset reticulate python environment, for more details, see:
  # https://github.com/greta-dev/greta/issues/444
  # Sys.unsetenv("RETICULATE_PYTHON")

  if (interactive()) configure_environment()
  else configure_environment(use_default=T)

}
