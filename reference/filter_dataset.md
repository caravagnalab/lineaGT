# Filters the input dataset.

Function used to filter observations, i.e. ISs, in the input dataframe
with coverage values.

## Usage

``` r
filter_dataset(
  cov.df,
  min_cov = 5,
  min_frac = 0.05,
  k_interval = c(10, 20),
  metric = "calinski_harabasz_score",
  seed = 5
)
```

## Arguments

- cov.df:

  Input coverage dataset. It must have at least the columns `coverage`,
  `timepoints`, `lineage`, `IS`, with the coverage values, timepoint,
  lineage and IS, respectively.

- min_cov:

  add

- min_frac:

  add

- k_interval:

  add

- metric:

  add

- seed:

  add

## Value

a dataset of the same shape as the input one, with filtered
observations.
