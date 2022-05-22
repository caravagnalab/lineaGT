#' Configure the reticulate environment
#'
#' @description Function to configure the Python dependencies on R.
#' If a Python environment is not available, the function will check if there is a version of
#' \code{conda} or \code{miniconda}, otherwise it will install \code{miniconda}, on which
#' install the Python package \code{pylineaGT}.
#'
#' @param env_name name of the \code{conda} environment to use, if available.
#'
#' @importFrom reticulate import conda_create conda_install install_miniconda miniconda_path
#' @export configure_environment

configure_environment = function(env_name="lineagt_env") {
  tryCatch(
    # try to import the python package
    expr = { reticulate::import("pylineaGT") },
    error = function(e) {
      tryCatch(
        # try to install the package in the already created environment
        expr = { reticulate::conda_install(env_name, "pylineaGT", pip=TRUE) },
        error = function(e) {
          # try to create the virtual environment
          tryCatch(expr = { reticulate::conda_create("lineagt_env") },
                   error = function (e) {
                     # install miniconda and the virtual environment
                     reticulate::install_miniconda()
                     reticulate::conda_create("lineagt_env") },
                   # install the python package
                   finally = reticulate::conda_install("lineagt_env", "pylineaGT", pip=TRUE)
          )
          }
      )
      }
    )

  reticulate::use_condaenv("lineagt_env", required=TRUE)
}
