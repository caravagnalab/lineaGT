#' Creates an object of class \code{mvnmm}.
#'
#' @description Function to fit the input data.
#'
#' @param dataset Input dataset. It has to be of shape \code{NxT}, where \code{N} is the number of
#' observations and \code{T} is the dimensionality of the dataset.
#' @param lineages List of the lineages into consideration.
#' @param columns List of columns to be selected from the dataset. If not specified,
#' all columns will be considered.
#' @param IS_values List of Insertion Sites identifier. It must be a list of length \code{N}, hence
#' for each observation a value is needed. If not specified, a general identifier \code{IS.1,...IS.N}
#' will be generated.
#' @param k_interval Interval of K values to test.
#' @param n_runs Number of runs to perform for each K.
#' @param steps Maximum number of steps for the inference.
#' @param lr Learning rate used in the inference.
#' @param covariance Covariance type for the Multivariate Gaussian.
#' @param random_state Value of the seed.
#' @return a \code{mvnmm} object, containing the input dataset, annotated with IS_values, N, K, T
#' specific of the dataset, the input IS and column names, a list params that will contain the
#' inferred parameters, the python object
#'
#' @importFrom dplyr filter mutate select group_by inner_join rename_with case_when all_of ungroup
#' @importFrom tidyr separate unite pivot_wider pivot_longer tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace_all
#' @importFrom MASS mvrnorm
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom reshape2 melt dcast
#' @importFrom grDevices gray
#' @importFrom stats coef lm runif
#' @importFrom utils combn read.csv
#' @importFrom VIBER variational_fit choose_clusters
#'
#' @export fit

fit = function(dataset, lineages, columns=list(), IS_values=list(), k_interval=c(10,30),
               n_runs=3, steps=500, lr=0.005, covariance="diag", random_state=25) {
  ic = data.frame(matrix(nrow=0, ncol=4)); colnames(ic) = c("BIC", "AIC", "ICL", "NLL")
  losses = data.frame(matrix(nrow=0, ncol=steps)); colnames(losses) = paste("iter_", 1:steps, sep="")
  grads = data.frame(matrix(nrow=0, ncol=steps)); colnames(grads) = paste("iter_", 1:steps, sep="")

  for (k in as.integer(k_interval[1]):as.integer(k_interval[2])) {
    for (run in 1:n_runs) {
      x_k = single_fit(k, dataset, lineages, columns, IS_values, run, steps, covariance, lr, random_state)

      kk = x_k$K
      n_iter = x_k$n_iter

      ic[paste(kk, run,sep=":"),] = x_k$IC
      losses[paste(kk, run, sep=":"),1:n_iter] = x_k$losses
      grads[paste(kk, run, "mean", sep=":"), 1:n_iter] = x_k$gradients$mean
      grads[paste(kk, run, "sigma", sep=":"), 1:n_iter] = x_k$gradients$sigma
      grads[paste(kk, run, "weights", sep=":"), 1:n_iter] = x_k$gradients$weights

      gc()
    }
  }
  selection = list("ic"=ic, "losses"=losses, "grads"=grads)

  best_k = get_best_k(selection, method="BIC")
  x = single_fit(best_k, dataset, lineages, columns, IS_values, steps=steps, lr=lr)
  x$runs = selection

  return(x)
}


single_fit = function(k, df, lineages, columns=list(), IS=list(), run=NULL, steps=500,
                      covariance="diag", lr=0.001, random_state=25) {
  print(paste("RUN", run, "- K =", k))
  x = initialize_object(k, df, lineages=lineages, columns=columns, IS_values=IS)
  x = run_inference(x, steps=as.integer(steps), covariance=covariance, lr=as.numeric(lr),
                    random_state=random_state)
  x = classifier(x)

  x$IC = compute_IC(x$py_model)

  x$losses = load_losses(x$py_model)
  x$gradients = load_params_gradients(x$py_model)

  x$n_iter = x$py_model$losses_grad_train$losses %>% length

  return(x)
}
