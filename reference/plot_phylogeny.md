# Clonal evolution trees

Clonal evolution trees

## Usage

``` r
plot_phylogeny(x, show_best = 1, min_frac = 0, highlight = c())
```

## Arguments

- x:

  a `mvnmm` object.

- show_best:

  the number of trees to visualize based on the computed score.

- min_frac:

  value in `[0,1]` representing the minimum abundance to show the
  clusters.

- highlight:

  a list of labels ID to show.

## Value

A list of `ggplot` objects with the estimated clone trees.
