# Input formats

``` r
library(lineaGT)
#> ✔ Loading ctree, 'Clone trees in cancer'. Support : <https://caravagn.github.io/ctree/>
#> ✔ Loading VIBER, 'Variational inference for multivariate Binomial mixtures'. Support : <https://caravagn.github.io/VIBER/>
#> ✔ Loading lineaGT, 'Lineage inference from gene therapy'. Support : <https://caravagnalab.github.io/lineaGT/>
# library(magrittr)
```

The `fit` function of the library requires two input datasets: one
dataset with the coverage observation and one with the mutations variant
allele number of reads and per-locus depth.

## Coverage Dataframe

The first dataframe requires the following columns:

- `IS`: the integration site ID,

- `timepoints`: the longitunal timepoint,

- `lineage`: the cell lineage name,

- `coverage` the number of reads assigned to the ISs.

If `lineage` and `timepoints` columns are not present, a single
longitunal observation and single lineage will be assumed.

A dataset example is the following:

``` r
data(cov.df.example)
cov.df.example
#> # A tibble: 428 × 4
#>    IS    timepoints lineage coverage
#>    <chr> <chr>      <chr>      <int>
#>  1 IS1   t1         l1           124
#>  2 IS1   t2         l1           190
#>  3 IS1   t1         l2             2
#>  4 IS1   t2         l2             6
#>  5 IS10  t1         l1             4
#>  6 IS10  t2         l1            14
#>  7 IS10  t1         l2             0
#>  8 IS10  t2         l2            12
#>  9 IS100 t1         l1             0
#> 10 IS100 t2         l1           418
#> # ℹ 418 more rows
```

## Mutations Dataframe

The second dataframe requires the following columns:

- `IS`: the integration site ID,

- `mutation`: the mutation ID,

- `timepoints`: the longitunal timepoint,

- `lineage`: the cell lineage name,

- `alt`: the per-locus variant allele number of reads,

- `dp`: the per-locus total number of reads, hence the per-locus depth.

A dataframe example is the following:

``` r
data(vaf.df.example)
vaf.df.example
#> # A tibble: 116 × 6
#>    IS    mutation timepoints lineage   alt    dp
#>    <chr> <chr>    <chr>      <chr>   <dbl> <dbl>
#>  1 IS1   mut1     t1         l1         58   134
#>  2 IS1   mut1     t2         l1         53   173
#>  3 IS1   mut1     t1         l2          0    26
#>  4 IS1   mut1     t2         l2         90   322
#>  5 IS10  mut12    t1         l1          0   372
#>  6 IS10  mut12    t2         l1          0   146
#>  7 IS10  mut12    t1         l2          0    57
#>  8 IS10  mut12    t2         l2          0   482
#>  9 IS11  mut13    t1         l1        160   372
#> 10 IS11  mut13    t2         l1          0   146
#> # ℹ 106 more rows
```
