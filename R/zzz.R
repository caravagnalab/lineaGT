.onLoad <- function(libname, pkgname) {
  Sys.unsetenv("RETICULATE_PYTHON")

  tryCatch(
    expr = reticulate::use_condaenv("lineagt-env", required = TRUE),
    error = function(e) NULL
  )

}
