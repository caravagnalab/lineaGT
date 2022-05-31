#' Creates an object of class \code{mvnmm}.
#'
#' @description Function to fit the input data.
#'
#' @param cov.df Input coverage dataset. It must have at least the columns \code{coverage}, \code{timepoints},
#' \code{lineage}, \code{IS}, with the coverage values, timepoint, lineage and IS, respectively.
#' @param vaf.df Input VAF dataset. It must have at least the columns \code{alt}, \code{dp}, \code{vaf},
#' \code{timepoints}, \code{lineage}, \code{IS}, \code{mutation}, with the number of reads for the mutated allele,
#' overall depth, vaf values, timepoint, lineage, IS and mutation, respectively.
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
#' @importFrom dplyr filter mutate select group_by inner_join rename_with case_when all_of ungroup slice
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

fit = function(cov.df,
               vaf.df=NULL,
               k_interval=c(10,30),
               n_runs=3,
               steps=500,
               lr=0.005,
               p=0.01,
               convergence=TRUE,
               covariance="diag",
               min_frac=0,
               show_progr=TRUE,
               store_grads=TRUE,
               store_losses=TRUE,
               random_state=25) {

  py_pkg = reticulate::import("pylineaGT")

  out = py_pkg$run_inference$Run(cov_df=cov.df %>% long_to_wide_cov(),
                                 lineages=cov.df$lineage %>% unique(),
                                 k_interval=list(as.integer(k_interval[1]), as.integer(k_interval[2])),
                                 n_runs=as.integer(n_runs),
                                 steps=as.integer(steps),
                                 lr=as.numeric(lr),
                                 p=as.numeric(p),
                                 convergence=convergence,
                                 covariance=covariance,
                                 show_progr=show_progr,
                                 store_grads=store_grads,
                                 store_losses=store_losses,
                                 random_state=as.integer(random_state))

  selection = list("ic"=out[[1]], "losses"=out[[2]], "grads"=out[[3]])

  best_k = get_best_k(selection, method="BIC")
  x = fit_singleK(best_k, cov.df, steps=steps, lr=lr, py_pkg)
  x$runs = selection

  if (!is.null(vaf.df)) x = run_viber(x, vaf.df, min_frac=min_frac)

  return(x)
}


