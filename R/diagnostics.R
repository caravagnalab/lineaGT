get_best_k = function(selection, method="BIC") {
  best = get_ic_df(selection) %>%
    group_by(K, variable) %>%
    dplyr::summarise(value_best=min(value))

  best_k = ((best %>% filter(variable=="BIC"))[(best %>% filter(variable=="BIC"))$value_best %>%
                                                 which.min(),"K"] %>%
              as.data.frame())[1,"K"] %>% droplevels() %>% levels() %>% as.integer()

  return(best_k)
}


get_ic_df = function(selection) {
  ic_df = selection$ic %>% rownames_to_column() %>%
    reshape2::melt(id=c("rowname")) %>%
    separate("rowname", into=c("K", "run")) %>%
    mutate(K=as.integer(K)) %>%
    mutate(K=factor(K, levels=paste(sort(K %>% unique()))))

  return(ic_df)
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
