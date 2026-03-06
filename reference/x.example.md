# Example mutation data

Example mutation data

## Usage

``` r
data(vaf.df.example)
```

## Format

An object of class `mvnmm` of length 20.

## Examples

``` r
data(vaf.df.example)
head(vaf.df.example)
#> # A tibble: 6 × 6
#>   IS    mutation timepoints lineage   alt    dp
#>   <chr> <chr>    <chr>      <chr>   <dbl> <dbl>
#> 1 IS1   mut1     t1         l1         58   134
#> 2 IS1   mut1     t2         l1         53   173
#> 3 IS1   mut1     t1         l2          0    26
#> 4 IS1   mut1     t2         l2         90   322
#> 5 IS10  mut12    t1         l1          0   372
#> 6 IS10  mut12    t2         l1          0   146
```
