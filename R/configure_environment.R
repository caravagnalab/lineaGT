#' Configure the reticulate environment
#'
#' @description Function to configure the Python dependencies on R.
#' If a Python environment is not available, the function will check if there is a version of
#' \code{conda} or \code{miniconda}, otherwise it will install \code{miniconda}, on which
#' install the Python package \code{pylineaGT}.
#'
#' @param env_name name of the \code{conda} environment to use, if available.
#'
#' @importFrom reticulate import conda_create conda_install install_miniconda miniconda_path conda_binary
#' @export configure_environment


configure_environment = function() {
  # install miniconda if needed
  if (!have_conda()) { install_miniconda_lineagt() }

  if (!have_lineagt_conda_env()) { create_conda_env_lineagt() }

  install_python_deps()
}


# find out whether the user has conda installed and visible
have_conda = function() {
  conda_bin = tryCatch(reticulate::conda_binary("auto"),
                       error = function(e) NULL)
  !is.null(conda_bin)
}


install_miniconda_lineagt = function() {
  reticulate::install_miniconda()
}


have_lineagt_conda_env = function(){
  tryCatch(expr = "lineagt-env" %in% reticulate::conda_list()$name,
           error = function(e) FALSE )
}


create_conda_env_lineagt = function() {
  reticulate::conda_create(envname="lineagt-env")
}


install_python_deps = function() {
  reticulate::conda_install(envname="lineagt-env", packages=c("pylineaGT"), pip=TRUE)
}


use_lineagt_conda_env = function() {
  tryCatch(
    expr = reticulate::use_condaenv("lineagt-env", required=TRUE),
    error = function(e) NULL
  )
}


using_lineagt_conda_env = function() {
  config = reticulate::py_discover_config()
  grepl("lineagt-env", config$python)
}



# configure_environment = function(env_name="lineagt-env") {
#   tryCatch(
#     # try to import the python package
#     expr = { reticulate::import("pylineaGT") },
#     error = function(e) {
#       tryCatch(
#         # try to install the package in the already created environment
#         expr = { reticulate::conda_install(env_name, "pylineaGT", pip=TRUE) },
#         error = function(e) {
#           # try to create the virtual environment
#           tryCatch(expr = { reticulate::conda_create("lineagt-env") },
#                    error = function (e) {
#                      # install miniconda and the virtual environment
#                      reticulate::install_miniconda()
#                      reticulate::conda_create("lineagt-env") },
#                    # install the python package
#                    finally = reticulate::conda_install("lineagt-env", "pylineaGT", pip=TRUE)
#           )
#           }
#       )
#       }
#     )
# }
