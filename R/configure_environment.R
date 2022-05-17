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
#'
#' @export configure_environment

configure_environment = function(env_name="r-reticulate") {
  tryCatch(expr = { reticulate::conda_install(env_name, "pylineaGT", pip=TRUE) },
           error = function(e) {
             tryCatch(expr = {
               reticulate::conda_create("r-reticulate")
               },
               error = function (e) { reticulate::install_miniconda(path=reticulate::miniconda_path(), update=TRUE, force=TRUE) },
               finally = reticulate::conda_install("r-reticulate", "pylineaGT", pip=TRUE)
             )
             }
           )
}
