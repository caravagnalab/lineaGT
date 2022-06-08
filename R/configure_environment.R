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


configure_environment = function(envname="lineagt-env", use_default=F) {
  # install miniconda if needed
  check_conda(use_default=use_default)

  envname = check_conda_env(envname=envname, use_default=use_default)

  check_python_deps(envname=envname, use_default=use_default)
}


check_conda = function(use_default=F) {
  if (!have_conda()) {
    cat("It is not prossibe to find a Anaconda or Miniconda version installed.\n")

    if (use_default) answ = "yes"
    else {
      cat("The Miniconda installer will be downloaded and used to install Miniconda. Proceed?\n")
      cat("(yes/no)\n")
      answ = readline()
    }

    if (answ == "yes") install_miniconda_lineagt()
    else { cat("Miniconda will not be installed."); return() }
  }
}


check_conda_env = function(envname="lineagt-env", use_default=F) {
  if (have_loaded_env()) {
    envname = sapply(reticulate::conda_list()$name, grepl, reticulate::py_discover_config()$python) %>%
      which() %>% names()
    return(envname)
  }

  if (!have_conda_env("lineagt-env")) {
    cat("The environment 'lineagt-env' is not present.\n")

    if (use_default) answ = "create"
    else {
      cat("Do you want to load an existing environment, to create a new one named 'lineagt-env' or to cancel?\n")
      cat("(load/create/cancel)\n")
      answ = readline()
      }

    if (answ == "create") {
      envname = "lineagt-env"
      create_conda_env()
    } else if (answ == "load") {
      cat("Insert the environment name: \n")
      envname = readline()
    } else {
      cat("No environment will be loaded nor created.")
      return()
    }
  } else {
    # cat("The environment 'lineagt-env' is already present and will be loaded!\n")
    envname = "lineagt-env"
  }

  use_conda_env(envname)

  return(envname)
}


check_python_deps = function(envname="lineagt-env", use_default=F) {
  try(install_python_deps(envname), silent=T)

  # if (!have_python_deps(envname)) {
  #   cat(paste("The required Python packages are not installed in the '", envname, "' environment.\n", sep=""))
  #
  #   if (use_default) answ = "yes"
  #   else {
  #     cat("Proceed with the installation?\n")
  #     cat("(yes/no)\n")
  #     answ = readline()
  #   }
  #
  #   if (answ == "yes") install_python_deps(envname)
  #   else return()
  # }
}

