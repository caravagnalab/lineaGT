.onLoad <- function(libname, pkgname) {
  tryCatch(expr = { reticulate::use_condaenv("lineagt_env", required=TRUE) },
           error = function(e) configure_environment()
           )
}
