# Fit the mutations clustering

Fit the mutations clustering

## Usage

``` r
fit_mutations(
  x,
  vaf.df = NULL,
  infer_phylo = TRUE,
  min_frac = 0,
  max_IS = NULL,
  highlight = list()
)
```

## Arguments

- x:

  a mvnmm object.

- vaf.df:

  dataframe with mutations data.

- infer_phylo:

  a Boolean indicating whether to infer also the phylogenetic evolution
  per cluster of ISs.

- min_frac:

  add

- max_IS:

  add

- highlight:

  add

## Value

A `mvnmm` object with the additional list `x.muts` containing the
estimated subclones from somatic mutations.
