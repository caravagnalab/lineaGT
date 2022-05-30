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


configure_environment = function(envname="lineagt-env") {
  if (interactive()) {
    # install miniconda if needed
    check_conda()

    envname = check_conda_env(envname=envname)

    check_python_deps(envname=envname)
  }

}


check_conda = function() {
  if (!have_conda()) {
    cat("It is not prossibe to find a Anaconda or Miniconda version installed.\n")
    cat("The Miniconda installer will be downloaded and used to install Miniconda. Proceed?\n")
    cat("(yes/no)\n")
    answ = readline()

    if (answ == "yes") install_miniconda_lineagt()
    else { cat("Miniconda will not be installed."); return() }
  }
}


check_conda_env = function(envname="lineagt-env") {
  if (!have_conda_env("lineagt-env")) {
    cat("The environment 'lineagt-env' is not present.\n")
    cat("Do you want to load an existing environment, to create a new one named 'lineagt-env' or to cancel?\n")
    cat("(load/create/cancel)\n")
    answ = readline()

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
    cat("The environment 'lineagt-env' is already present!\n")
    cat("Do you want to load it or do you want to load another environment 'envname'?\n")
    cat("(load/envname)\n")
    answ = readline()
    if (answ != "load") envname = readline()
  }

  use_conda_env(envname)

  return(envname)
}


check_python_deps = function(envname="lineagt-env") {
  if (!have_python_deps(envname)) {
    cat(paste("The required Python packages are not installed in the '", envname, "' environment.\n", sep=""))
    cat("Proceed with the installation?\n")
    cat("(yes/no)\n")
    answ = readline()

    if (answ == "yes") install_python_deps(envname)
    else return()
  }
}


have_python_deps = function(envname="lineagt-env", packages="pylineagt") {
    tryCatch(expr = packages %in% reticulate::py_list_packages(envname)$package,
             error = function(e) FALSE)
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


have_conda_env = function(envname="lineagt-env"){
  tryCatch(expr = envname %in% reticulate::conda_list()$name,
           error = function(e) FALSE )
}


create_conda_env = function(envname="lineagt-env") {
  reticulate::conda_create(envname=envname)
}


install_python_deps = function(envname="lineagt-env") {
  reticulate::conda_install(envname=envname, packages=c("pylineaGT"), pip=TRUE)
}


use_conda_env = function(envname="lineagt-env") {
  Sys.unsetenv("RETICULATE_PYTHON")
  tryCatch(
    expr = reticulate::use_condaenv(envname, required=TRUE),
    error = function(e) NULL
  )
}


using_conda_env = function(envname="lineagt-env") {
  config = reticulate::py_discover_config()
  grepl(envname, config$python)
}

