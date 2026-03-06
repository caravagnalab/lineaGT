# Extract the estimated posterior probabilities.

Returns a dataframe of shape `NxK` with the posterior distribution
`p(k|n)` for each observation `n` to belong to cluster `k`.

## Usage

``` r
get_z_probs(x)
```

## Arguments

- x:

  a mvnmm object.

## Value

dataframe of posterior distributions.

## Examples

``` r
if (FALSE) get_z_probs(x)
```
