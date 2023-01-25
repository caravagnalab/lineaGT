#' Creates an object of class \code{mvnmm}.
#'
#' @description Function to fit the input data.
#'
#' @param cov.df Input coverage dataset. It must have at least the columns \code{coverage}, and \code{IS},
#' and additional columns \code{timepoints} and \code{lineage}, will be added if missing, assuming single
#' timepoint and lineage.
#' @param vaf.df Input VAF dataset. If not \code{NULL}, the mutations clustering will be performed.
#' It must have at least the columns \code{mutation}, \code{IS}, \code{alt}, \code{dp}, and additional
#' \code{vaf}, \code{timepoints}, \code{lineage}, \code{IS}, \code{mutation}, with the number of reads for the mutated allele,
#' overall depth, vaf values, timepoint, lineage, IS and mutation, respectively.
#' @param infer_phylogenies A Boolean. If set to \code{TRUE}, the function will also compute and attach to
#' the returned object the phylogenetic trees for each cluster.
#' @param k_interval Interval of K values to test.
#' @param n_runs Number of runs to perform for each K.
#' @param steps Maximum number of steps for the inference.
#' @param lr Learning rate used in the inference.
#' @param p Numeric value used to check the convergence of the parameters.
#' @param min_frac add
#' @param max_IS add
#' @param check_conv A Boolean. If set to \code{TRUE}, the function will check for early convergence,
#' otherwise it will perform \code{steps} iterations.
#' @param covariance Covariance type for the Multivariate Gaussian.
#' @param hyperparams add
#' @param timepoints_to_int add
#' @param show_progr A Boolean. If \code{TRUE}, the progression bar will be shown during inference.
#' @param store_grads A Booolean. If \code{TRUE}, the gradient norms for the parameters at each
#' iteration will be stored.
#' @param store_losses A Boolean. If \code{TRUE}, the computed losses for the parameters at each
#' iteration will be stored.
#' @param store_params A Boolean. If \code{TRUE}, the estimated parameters at each iteration will be stored.
#' @param seed_optim add
#' @param seed Value of the seed.
#' @param sample_id add
#' @param default_lm add
#'
#' @return A \code{mvnmm} object, containing the input dataset, annotated with IS_values, N, K, T
#' specific of the dataset, the input IS and column names, a list params that will contain the
#' inferred parameters, the python object
#'
#' @importFrom dplyr filter mutate select group_by inner_join rename_with
#' @importFrom dplyr case_when all_of ungroup slice pull
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
#' @importFrom ctree ctrees
#'
#' @export fit

# hyperparams A list of hyperparams to set. Available values are:
# \code{mean_loc}, the center of the mean prior (default set to the sample mean),
# \code{mean_scale}, the variance of the mean prior (default set to the sample variance),
# \code{var_scale}, the variance of the variance prior (default set to \code{400}),
# \code{max_var}, the maximum value attributed to the variance (no default value),
# \code{min_var} the minimum value attributed to the variance (no default limit).

fit = function(cov.df,
               vaf.df=NULL,
               infer_phylogenies=TRUE,
               infer_growth=TRUE,
               k_interval=c(5,15),
               n_runs=1,
               steps=500,
               min_steps=20,
               lr=0.005,
               p=1,
               min_frac=0,
               max_IS=NULL,
               check_conv=TRUE,
               covariance="full",
               hyperparams=list(),
               default_lm=TRUE,
               timepoints_to_int=list(),
               show_progr=FALSE,
               store_grads=TRUE,
               store_losses=TRUE,
               store_params=FALSE,
               seed_optim=TRUE,
               seed=6,
               seed_init=reticulate::py_none(),
               sample_id="") {

  py_pkg = reticulate::import("pylineaGT")

  cov.df = cov.df %>% check_cov_dimensions()

  max_k = cov.df %>% check_max_k()
  k_interval = check_k_interval(k_interval, max_k)

  cli::cli_process_start("\nStarting lineaGT model selection to retrieve the optimal number of clones\n")
  out = py_pkg$run$run_inference(cov_df=cov.df %>% long_to_wide_cov(),
                                 lineages=cov.df$lineage %>% unique(),
                                 k_interval=list(as.integer(k_interval[1]), as.integer(k_interval[2])),
                                 n_runs=as.integer(n_runs),
                                 steps=as.integer(steps),
                                 lr=as.numeric(lr),

                                 check_conv=check_conv,
                                 min_steps=as.integer(min_steps),
                                 p=as.numeric(p),

                                 default_lm=default_lm,
                                 covariance=covariance,
                                 hyperparams=reticulate::py_dict(keys=names(hyperparams),
                                                                     values=as.numeric(hyperparams)),
                                 show_progr=show_progr,
                                 store_grads=store_grads,
                                 store_losses=store_losses,
                                 store_params=store_params,

                                 seed_optim=seed_optim,
                                 seed=as.integer(seed)
                                 )
  cli::cli_process_done()

  selection = list("ic"=out[[1]], "losses"=out[[2]], "grads"=out[[3]], "params"=out[[4]]) %>%
    get_selection_df()

  best_k = get_best_k(selection, method="BIC")$K
  best_init_seed = get_best_k(selection, method="BIC")$init_seed
  best_seed = get_best_k(selection, method="BIC")$seed

  cli::cli_process_start("Fitting model to cluster ISs")
  x = fit_singleK(k=best_k,
                  cov.df=cov.df,
                  steps=steps,

                  hyperparams=hyperparams,
                  default_lm=default_lm,
                  covariance=covariance,

                  lr=lr,
                  p=p,
                  min_steps=min_steps,
                  check_conv=check_conv,

                  store_params=store_params,

                  seed_optim=FALSE,
                  seed=best_seed,
                  init_seed=best_init_seed,

                  timepoints_to_int=timepoints_to_int,
                  py_pkg=py_pkg)

  cli::cli_process_done(msg_done=paste0("Found ", x$K, " clones of ISs!"))

  x$runs = selection

  if (!is.null(vaf.df)) {
    cli::cli_process_start("Fitting model to cluster mutations")
    x = fit_mutations(x, vaf.df, infer_phylo=infer_phylogenies, min_frac=min_frac, max_IS=max_IS)
    cli::cli_process_done()
  }

  if (have_muts_fit(x)) mutations = T else mutations = F
  x$population.df = get_muller_pop(x, mutations=mutations, timepoints_to_int=timepoints_to_int, estimate_npops=F)

  if (infer_growth) {
    cli::cli_process_start("Fitting model to estimate population growth rates")
    x = fit_growth_rates(x, steps=steps, timepoints_to_int=timepoints_to_int, py_pkg=py_pkg)
    cli::cli_process_done()
  }

  x$sample_id = sample_id

  return(x)
}


check_cov_dimensions = function(cov.df) {
  tp = cov.df %>% dplyr::pull(timepoints) %>% unique() %>% as.character()
  lin = cov.df %>% dplyr::pull(lineage) %>% unique() %>% as.character()

  # all the possible combinations of timepoints.lineage
  dimensions = apply(expand.grid(tp, lin), 1, paste, collapse=".")

  # find the combinations of timepoints.lineage present in the data
  combs = cov.df %>%
    long_to_wide_cov() %>%
    dplyr::select(starts_with("cov")) %>%
    colnames() %>%
    str_replace_all("cov.","")

  # find the combinations not present in the data
  missing.vals = setdiff(dimensions, combs) %>%
    tibble::as_tibble() %>%
    tidyr::separate(value, into=c("timepoints", "lineage"), sep="[.]") %>%
    dplyr::mutate(coverage=0) %>%
    dplyr::mutate(IS=cov.df$IS[1])

  if (nrow(missing.vals) == 0)
    return(cov.df %>% dplyr::ungroup())

  try(expr = {
    missing.vals = missing.vals %>%
      dplyr::mutate(timepoints=as.integer(timepoints))
  }, silent=T)

  return(
    cov.df %>%
      dplyr::ungroup() %>%
      dplyr::add_row(
        missing.vals
      )
    )

}


check_max_k = function(cov.df) {
  # function to check the maximum value of K
  # it has to be set to the number of distinct observations -1
  # otherwise the Kmeans will give an error
  n_dinstinct = cov.df %>%
    long_to_wide_cov() %>%
    dplyr::select(dplyr::starts_with("cov")) %>%
    unique() %>%
    nrow()

  return(
    n_dinstinct - 1
  )
}


check_k_interval = function(k_interval, max_k) {
  if (k_interval[2] > max_k) {
    k_interval[2] = max_k
    k_interval[1] = max(2, k_interval[2] - 10)
  }

  return(k_interval)
}

