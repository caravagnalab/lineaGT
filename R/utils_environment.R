have_loaded_env = function() {
  return(grepl("envs", reticulate::py_discover_config()$python))
}


have_conda_env = function(envname="lineagt-env"){
  tryCatch(expr = envname %in% reticulate::conda_list()$name,
           error = function(e) FALSE )
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


create_conda_env = function(envname="lineagt-env") {
  reticulate::conda_create(envname=envname)
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


install_miniconda_lineagt = function() {
  reticulate::install_miniconda()
}


install_python_deps = function(envname="lineagt-env") {
  reticulate::conda_install(envname=envname, packages=c("pylineaGT"), pip=TRUE)
}
