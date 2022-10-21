.onLoad <- function(libname, pkgname) {

  pk = 'lineaGT'
  pk_l = 'Lineage inference from gene therapy'
  www = "https://caravagnalab.github.io/lineaGT/"

  cli::cli_alert_success(
    'Loading {.field {pk}}, {.emph \'{pk_l}\'}. Support : {.url { www}}' )

  # unset reticulate python environment, for more details, see:
  # https://github.com/greta-dev/greta/issues/444
  # Sys.unsetenv("RETICULATE_PYTHON")

  # configure_environment()

  if (interactive()) configure_environment(use_default=F)
  else configure_environment(use_default=T)
}
