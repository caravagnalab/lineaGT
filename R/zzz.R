.onLoad <- function(libname, pkgname) {
  tryCatch(expr = { reticulate::import("pylineaGT") },
           error = function(e) {
             tryCatch(expr = { reticulate::conda_install("r-reticulate", "pylineaGT", pip=TRUE) },
                      error = function(e) {
                        reticulate::install_miniconda(path=reticulate::miniconda_path(), update=TRUE, force=TRUE)
                        reticulate::conda_install("r-reticulate", "pylineaGT", pip=TRUE)
                      })
             } )
}
