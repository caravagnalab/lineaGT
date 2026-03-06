# Visualize the regression given the infered growth rates.

Function to visualize the growht of each clone and lineage.

## Usage

``` r
plot_growth_regression(
  x,
  highlight = c(),
  min_frac = 0,
  mutations = F,
  timepoints_to_int = list(),
  fit = F,
  show_best = T,
  ratio = NULL
)
```

## Arguments

- x:

  a mvnmm object.

- highlight:

  a vector of clusters IDs to highlight in the plot.

- min_frac:

  min_frac numeric value in `[0,1]` representing the minimum abundance
  to highlight a clone.

- mutations:

  Boolean. If set to `TRUE`, the growth will be visualize for each
  cluster of mutations.

- timepoints_to_int:

  a list to map each `timepoint` value to an integer.

- fit:

  add

## Examples

``` r
if (FALSE) plot_exp_fit(x)
```
