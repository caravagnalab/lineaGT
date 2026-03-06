# Extract the estimated Cholesky matrices, used to factorise the covariance matrix.

Returns a list with `K` dataframes, each of dimension `TxT`,
corresponding to the covariance matrices estimated for each clone `k`

## Usage

``` r
get_covariance_Cholesky(x)
```

## Arguments

- x:

  a mvnmm object.

## Value

list of the estimated covariance matrices.

## Examples

``` r
if (FALSE) get_covariance_Cholesky(x)
```
