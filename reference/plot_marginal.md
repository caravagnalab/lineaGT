# Histogram of the marginal distribution of each dimension

Function to plot the marginal distribution of the coverage, for each
timepoint, colored by cluster.

## Usage

``` r
plot_marginal(
  x,
  min_frac = 0,
  highlight = c(),
  binwidth = 10,
  show_dens = T,
  timepoints_to_int = list(),
  facet_lin = F,
  single_plot = F
)
```

## Arguments

- x:

  a mvnmm object.

- highlight:

  a vector of clusters IDs to highlight in the plot.

- binwidth:

  numeric value representing the histogram binwidth.

## Examples

``` r
if (FALSE) plot_marginal(x)
```
