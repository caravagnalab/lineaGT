# Infer growth rates for each clone and subclone.

Infer growth rates for each clone and subclone.

## Usage

``` r
fit_growth_rates(
  x,
  steps = 500,
  warmup_steps = 500,
  highlight = c(),
  timepoints_to_int = c(),
  timepoints = get_timepoints(x),
  growth_model = "exp.log",
  which = "pop",
  force = T,
  tree_score = 1,
  py_pkg = NULL,
  mutations = F
)
```

## Arguments

- x:

  a mvnmm object.

- steps:

  maximum number of steps for inference.

- highlight:

  set of clusters to run the inference for. If not specified, it will be
  run on all the clusters.

- timepoints_to_int:

  if the provided timepoints are not integers nor a timepoints-to-int
  list is stored in `x`, a list mapping their values to integers is
  required.

- growth_model:

  string specifying the type of growth model to use, between `exp` and
  `log` corresponding to exponential and logistic models, respectively.

- force:

  if the model has already been fitted, setting `force` to `FALSE` will
  keep the computed rates. Setting `force` to `TRUE` will fit the model
  again for the specified clusters.

## Value

A `mvnmm` object with the additional tibble `growth.rates` containing
the estimated population genetics parameters.
