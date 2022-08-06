get_selection_df = function(selection) {
  runs = list()

  runs$IC = get_IC(selection)
  runs$losses = get_losses(selection=selection, train=T)
  runs$grads = get_gradient_norms(selection)
  runs$params = get_params_train(selection$params)

  return(runs)

}


get_best_k = function(selection, method="BIC") {
  return(
    selection$IC %>%

      dplyr::mutate(K=as.integer(K), run=as.integer(run), value=as.numeric(value)) %>%
      dplyr::rename(mm=method) %>%
      dplyr::filter(mm==method) %>%
      dplyr::filter(value==min(value)) %>%

      tidyr::separate("id", into=c("K.init", "else"), sep="[.]") %>%
      dplyr::mutate(K.init=as.integer(K.init)) %>%
      dplyr::select(K, seed, init_seed) %>% unique() %>% as.list()
    )
}


get_IC = function(selection) {
  return(selection$ic %>%
           reshape2::melt(id=c("K","run","id","seed","init_seed"), variable.name="method") %>%
           tibble::as_tibble() %>%
           dplyr::mutate(K=as.integer(K), run=as.integer(run))
         )
}


get_losses = function(x=NULL, selection=NULL, train=FALSE) {
  if (train) return(selection$losses %>%
                     tibble::as_tibble() %>%
                     dplyr::mutate(K=as.integer(K), run=as.integer(run)))
  return(
    x$losses %>%
      tibble::as_tibble() %>%
      dplyr::rename(losses=value) %>%
      tibble::rownames_to_column(var="index") %>%
      dplyr::mutate(index=as.numeric(index))
           )
}


get_gradient_norms = function(selection) {
  return(selection$grads %>%
           tibble::as_tibble() %>%
           dplyr::mutate(K=as.integer(K),
                  run=as.integer(run),
                  param=str_replace_all(param, "_param","")))
}


get_params_train = function(params) {
  par.list = list()
  K = params$K %>% unlist()
  run = params$run %>% unlist()
  id = params$id %>% unlist()
  param.name = params$param
  value = params$params_values

  df = tibble::tibble("K"=K,
                      "run"=run,
                      "id"=id,
                      "param"=param.name,
                      "value"=lapply(value, tibble::as_tibble_col) %>% lapply(rownames_to_column, var="step"))

  return(df)
}


load_losses = function(py_model) {
  return(py_model$losses_grad_train$losses)
}


load_params_gradients = function(py_model) {
  gradients = list()
  gradients$mean = py_model$losses_grad_train$gradients$mean_param
  gradients$sigma = py_model$losses_grad_train$gradients$sigma_vector_param
  gradients$weights = py_model$losses_grad_train$gradients$weights_param
  return(gradients)
}


compute_IC = function(py_model) {
  IC = list()
  IC$BIC = py_model$compute_ic(method="BIC")$numpy()
  IC$AIC = py_model$compute_ic(method="AIC")$numpy()
  IC$ICL = py_model$compute_ic(method="ICL")$numpy()
  IC$NLL = py_model$nll$numpy()
  return(IC)
}
