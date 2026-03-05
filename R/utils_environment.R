have_loaded_env = function() {
  return(grepl("envs", reticulate::py_discover_config()$python))
}


have_python_deps = function(envname="", py_pkgs=c("pylineagt")) {
  if (envname == "")
    envname = which_conda_env()

  tryCatch(
    expr = {
      ll = py_pkgs %in% reticulate::py_list_packages(envname)$package
      names(ll) = py_pkgs
      return(ll)
    },
    error = function(e) FALSE)
}


which_conda_env = function() {
  if (have_loaded_env())
    return(
      sapply(reticulate::conda_list()$name, grepl, reticulate::py_discover_config()$python) %>%
        which() %>%
        names()
    )

  return(cat("No loaded environments!"))
}


load_conda_env = function(envname="lineagt-env") {
  Sys.unsetenv("RETICULATE_PYTHON")
  tryCatch(
    expr = reticulate::use_condaenv(envname, required=TRUE),
    error = function(e) {
      message("To change the loaded environment, you need to restart the R session!")
    }
  )
}


have_conda_env = function(envname="lineagt-env"){
  tryCatch(expr = envname %in% reticulate::conda_list()$name,
           error = function(e) FALSE )
}


# find out whether the user has conda installed and visible
have_conda = function() {
  conda_bin = tryCatch(reticulate::conda_binary("auto"),
                       error = function(e) NULL)
  !is.null(conda_bin)
}


create_conda_env = function(envname="lineagt-env") {
  reticulate::conda_create(envname=envname)
}


using_conda_env = function(envname="lineagt-env") {
  config = reticulate::py_discover_config()
  grepl(envname, config$python)
}


install_miniconda_lineagt = function() {
  reticulate::install_miniconda()
}


install_python_deps = function(envname="lineagt-env") {
  reticulate::conda_install(envname=envname, packages=c("pylineaGT"),
                            pip=TRUE,
                            pip_options=c("--timeout 200"))
}
