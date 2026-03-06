# 2D scatterplot and density

Function to plot the scatterplots of the coverage, one timepoint against
the other, together with the Gaussian density.

## Usage

``` r
plot_scatter_density(x, plot_density = T, highlight = c(), min_frac = 0)
```

## Arguments

- x:

  a mvnmm object.

- plot_density:

  a Boolean. If set to FALSE, the Gaussian density will not be
  displayed.

- highlight:

  a vector of clusters IDs to highlight in the plot.

- facet:

  a Boolean. If set to TRUE, the plot will be faceted on the labels.

## Value

a list of plots, one for each timepoints combination.

## Examples

``` r
if (FALSE) plots = plot_scatter_density(x)
```
