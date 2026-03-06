# Creates an object of class `mvnmm`.

Function to fit the input data.

## Usage

``` r
fit(
  cov.df,
  vaf.df = NULL,
  infer_phylogenies = TRUE,
  infer_growth = TRUE,
  k_interval = c(5, 15),
  n_runs = 1,
  steps = 500,
  min_steps = 20,
  lr = 0.005,
  p = 1,
  min_frac = 0,
  max_IS = NULL,
  check_conv = TRUE,
  covariance = "full",
  hyperparams = list(),
  default_lm = TRUE,
  timepoints_to_int = list(),
  show_progr = FALSE,
  store_grads = TRUE,
  store_losses = TRUE,
  store_params = FALSE,
  seed_optim = TRUE,
  seed = 6,
  seed_init = reticulate::py_none(),
  py_pkg = NULL,
  sample_id = ""
)
```

## Arguments

- cov.df:

  Input coverage dataset. It must have at least the columns `coverage`,
  and `IS`, and additional columns `timepoints` and `lineage`, will be
  added if missing, assuming single timepoint and lineage.

- vaf.df:

  Input VAF dataset. If not `NULL`, the mutations clustering will be
  performed. It must have at least the columns `mutation`, `IS`, `alt`,
  `dp`, and additional `vaf`, `timepoints`, `lineage`, `IS`, `mutation`,
  with the number of reads for the mutated allele, overall depth, vaf
  values, timepoint, lineage, IS and mutation, respectively.

- infer_phylogenies:

  A Boolean. If set to `TRUE`, the function will also compute and attach
  to the returned object the phylogenetic trees for each cluster.

- k_interval:

  Interval of K values to test.

- n_runs:

  Number of runs to perform for each K.

- steps:

  Maximum number of steps for the inference.

- lr:

  Learning rate used in the inference.

- p:

  Numeric value used to check the convergence of the parameters.

- min_frac:

  add

- max_IS:

  add

- check_conv:

  A Boolean. If set to `TRUE`, the function will check for early
  convergence, otherwise it will perform `steps` iterations.

- covariance:

  Covariance type for the Multivariate Gaussian.

- hyperparams:

  add

- default_lm:

  add

- timepoints_to_int:

  add

- show_progr:

  A Boolean. If `TRUE`, the progression bar will be shown during
  inference.

- store_grads:

  A Booolean. If `TRUE`, the gradient norms for the parameters at each
  iteration will be stored.

- store_losses:

  A Boolean. If `TRUE`, the computed losses for the parameters at each
  iteration will be stored.

- store_params:

  A Boolean. If `TRUE`, the estimated parameters at each iteration will
  be stored.

- seed_optim:

  add

- seed:

  Value of the seed.

- sample_id:

  add

## Value

A `mvnmm` object, containing the input dataset, annotated with
IS_values, N, K, T specific of the dataset, the input IS and column
names, a list params that will contain the inferred parameters, the
python object
