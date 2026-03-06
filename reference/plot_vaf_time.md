# VAF over time

Function to plot the VAFs of the mutations along time

## Usage

``` r
plot_vaf_time(x, min_frac = 0, highlight = c(), timepoints_to_int = list())
```

## Arguments

- x:

  a `mvnmm` object.

- min_frac:

  value in `[0,1]` representing the minimum abundance to show the
  clusters.

- highlight:

  a list of labels ID to show.

- timepoints_to_int:

  a list to map each `timepoint` value to an integer.

## Value

a `ggplot` object.

## Examples

``` r
if (FALSE) plot_vaf_time(x, min_frac=0.1)
```
