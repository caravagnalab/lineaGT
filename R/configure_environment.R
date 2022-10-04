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
      cli::cli_alert_warning("The Miniconda installer will be downloaded and used to install Miniconda. Proceed? (yes/no)\n")
      answ = readline()
    }

    if (answ == "yes") install_miniconda_lineagt()
    else { cli::cli_alert_info("Miniconda will not be installed."); return() }
  }
}


check_conda_env = function(envname="lineagt-env", use_default=F) {
  if (have_loaded_env()) {
    envname = sapply(reticulate::conda_list()$name, grepl, reticulate::py_discover_config()$python) %>%
      which() %>% names()
    cli::cli_alert_warning(paste0("The '", envname, "' environment is already loaded!"))
    return(envname)
  }

  if (!have_conda_env("lineagt-env")) {
    cli::cli_alert_info("The environment 'lineagt-env' is not present.\n")

    if (use_default) answ = "create"
    else {
      cli::cli_alert_warning("Do you want to load an existing environment, to create a new one named 'lineagt-env' or to cancel? (load/create/cancel)\n")
      answ = readline()
      }

    if (answ == "create") {
      envname = "lineagt-env"
      create_conda_env()
    } else if (answ == "load") {
      cli::cli_alert_info("Insert the environment name: ")
      envname = readline()
    } else {
      cli::cli_alert_info("No environment will be loaded nor created.")
      return()
    }

  } else {
    cli::cli_alert_info("The environment 'lineagt-env' is already present and will be loaded!\n")
    envname = "lineagt-env"
  }

  load_conda_env(envname)

  return(envname)
}


check_python_deps = function(envname="lineagt-env", use_default=F) {
  install_python_deps(envname)

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

