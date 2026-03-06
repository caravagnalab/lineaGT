# VAF 2D scatterplot

Function to plot the VAFs of the mutations one timepoint against the
other

## Usage

``` r
plot_vaf(x, min_frac = 0, highlight = c(), wrap = T)
```

## Arguments

- x:

  a `mvnmm` object.

- min_frac:

  value in `[0,1]` representing the minimum abundance to show the
  clusters.

- highlight:

  a list of labels ID to show.

- wrap:

  a Boolean to wrap the scatterplots of multiple lineages in a unique
  plot.

## Value

list of VAF scatterplots.

## Examples

``` r
if (FALSE) plot_vaf(x, min_frac=0.1)
```
