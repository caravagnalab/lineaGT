get_best_k = function(selection, method="BIC") {
  return(
    selection$ic %>%
      reshape2::melt(id=c("K","run","seed"), variable.name="mm") %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(K=as.integer(K), run=as.integer(run), value=as.numeric(value)) %>%
      dplyr::filter(mm==method) %>%
      # dplyr::group_by(mm, K) %>%
      # dplyr::summarise(mean_value=mean(value)) %>%
      # dplyr::ungroup() %>%
      # dplyr::filter(mean_value==min(mean_value)) %>%
      dplyr::filter(value==min(value)) %>%
      dplyr::select(K,seed) %>% as.list()
    )
}


get_IC = function(runs) {
  return(runs$ic %>%
           reshape2::melt(id=c("K","run","id","seed"), variable.name="method") %>%
           dplyr::as_tibble() %>%
           dplyr::mutate(K=as.integer(K), run=as.integer(run))
         )
}


get_losses = function(x=NULL, runs=NULL, train=FALSE) {
  if (train) return(runs$losses %>%
                     dplyr::as_tibble() %>%
                     dplyr::mutate(K=as.integer(K), run=as.integer(run)))
  return(
    x$losses %>%
      tibble::as_tibble() %>%
      dplyr::rename(losses=value) %>%
      tibble::rownames_to_column(var="index") %>%
      mutate(index=as.numeric(index))
           )
}


get_gradient_norms = function(runs) {
  return(runs$grads %>%
           dplyr::as_tibble() %>%
           dplyr::mutate(K=as.integer(K),
                  run=as.integer(run),
                  param=str_replace_all(param, "_param","")))
}


compute_IC = function(py_model) {
  IC = list()
  IC$BIC = py_model$compute_ic(method="BIC")$numpy()
  IC$AIC = py_model$compute_ic(method="AIC")$numpy()
  IC$ICL = py_model$compute_ic(method="ICL")$numpy()
  IC$NLL = py_model$nll$numpy()
  return(IC)
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


get_selection_df = function(selection) {
  runs = list()

  runs$IC = get_IC(selection)
  runs$losses = get_losses(runs=selection, train=T)
  runs$grads = get_gradient_norms(selection)
  runs$params = get_params_train(selection$params)

  return(runs)

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
