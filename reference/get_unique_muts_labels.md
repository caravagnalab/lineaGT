# Retrieve the list of unique labels of mutation clusters.

Function to retrieve the list of unique labels of mutations clusters, of
the form `C_c1.Cm1`, where `c1` is the clone identifier and `m1` is the
subclone identifier.

## Usage

``` r
get_unique_muts_labels(x, clusters = c(), verbose = F)
```

## Arguments

- x:

  a mvnmm object.

- clusters:

  a vector-like variable, with the identifiers of the clones we want to
  retrieve the subclone labels from. If empty, all the labels will be
  returned.

## Value

vector of mutations labels.

## Examples

``` r
if(FALSE) get_unique_muts_labels(x, c("C_0"))
```
