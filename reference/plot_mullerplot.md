# Muller plot

Function to visualize the mullerplot for the fitted object.

## Usage

``` r
plot_mullerplot(
  x,
  which = "frac",
  highlight = c(),
  min_frac = 0,
  min_abundance = 0,
  estimate_npops = FALSE,
  vcn = NULL,
  rm_mixt = FALSE,
  timepoints_to_int = c(),
  mutations = F,
  single_clone = T,
  tree_score = 1,
  legend.pos = "right"
)
```

## Arguments

- x:

  a mvnmm object.

- which:

  string among `"frac","pop","fitness"` determining whether to plot the
  coverage normalized in `[0,1]`, as absolute clone abundance, or with
  each clone colored by the growth rate, computed assuming a exponential
  growth.

- highlight:

  a vector of clusters IDs to highlight in the plot.

- min_frac:

  min_frac numeric value in `[0,1]` representing the minimum abundance
  to highlight a clone.

- rm_mixt:

  remove clusters estimated as polyclonal

- timepoints_to_int:

  a list to map each `timepoint` value to an integer.

- mutations:

  Boolean. If set to `TRUE`, also the clusters of mutations will be
  visualized.

- single_clone:

  Boolean. If `mutations` and `single_clone` are set to `TRUE`, only the
  clones reported in `highlight` and the respective subclones will be
  visualised.

- tree_score:

  add

- legend.pos:

  add

## Examples

``` r
if (FALSE) plot_mullerplot(x, wrap=T)
```
