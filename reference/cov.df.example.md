# Example coverage data

Example coverage data

## Usage

``` r
data(cov.df.example)
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 428
rows and 4 columns.

## Examples

``` r
data(cov.df.example)
head(cov.df.example)
#> # A tibble: 6 × 4
#>   IS    timepoints lineage coverage
#>   <chr> <chr>      <chr>      <int>
#> 1 IS1   t1         l1           124
#> 2 IS1   t2         l1           190
#> 3 IS1   t1         l2             2
#> 4 IS1   t2         l2             6
#> 5 IS10  t1         l1             4
#> 6 IS10  t2         l1            14
```
