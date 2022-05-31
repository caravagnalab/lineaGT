# get_best_k = function(selection, method="BIC") {
#   best = get_ic(selection) %>%
#     group_by(K, variable) %>%
#     dplyr::summarise(value_best=min(value))
#
#   best_k = ((best %>% filter(variable=="BIC"))[(best %>% filter(variable=="BIC"))$value_best %>%
#                                                  which.min(),"K"] %>%
#               as.data.frame())[1,"K"] %>% droplevels() %>% levels() %>% as.integer()
#
#   return(best_k)
# }


get_best_k = function(selection, method="BIC") {
  return(selection$ic %>%
           reshape2::melt(id=c("K","run"), variable.name="mm") %>%
           dplyr::as_tibble() %>%
           dplyr::mutate(K=as.integer(K), run=as.integer(run)) %>%
           dplyr::filter(mm==method) %>%
           dplyr::slice(which.min(value)) %>%
           .$K)
}


get_ic = function(selection) {
  return(selection$ic %>%
           reshape2::melt(id=c("K","run"), variable.name="method") %>%
           dplyr::as_tibble() %>%
           dplyr::mutate(K=as.integer(K), run=as.integer(run)))
  # ic_df = selection$ic %>% tibble::rownames_to_column(var="k:run") %>%
  #   reshape2::melt(id=c("k:run")) %>%
  #   separate("k:run", into=c("K", "run"), sep=":") %>%
  #   mutate(K=as.integer(K)) %>%
  #   mutate(K=factor(K, levels=paste(sort(K %>% unique()))))
  #
  # return(ic_df)
}


get_losses = function(selection) {
  return(selection$losses %>%
           dplyr::as_tibble() %>%
           dplyr::mutate(K=as.integer(K), run=as.integer(run)))
}


get_gradient_norms = function(selection) {
  return(selection$grads %>%
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
