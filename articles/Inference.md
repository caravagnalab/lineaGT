# Lineage inference

`lineaGT` Bayesian models are implemented in Python using the Pyro
probabilistic programming language and thus require a working Python
programme and `pylineaGT` package. We recommend reading the [Get
Started](https://caravagnalab.github.io/lineaGT/articles/lineaGT.html)
article to correctly install the Python package.

``` r
reticulate::conda_create(envname="lineaGT", python_version="3.8")
#> + /usr/share/miniconda/bin/conda create --yes --name lineaGT 'python=3.8' --quiet -c conda-forge
#> [1] "/usr/share/miniconda/envs/lineaGT/bin/python"
reticulate::use_condaenv("lineaGT")
reticulate::conda_install(envname="lineaGT", packages="pylineaGT", pip=TRUE)
#> [1] "pylineaGT"
py_pkg = reticulate::import("pylineaGT")
```

Remember to load the new environment in order to have the `pylineaGT`
module accessible.

``` r
library(lineaGT)
#> Warning: replacing previous import 'cli::num_ansi_colors' by
#> 'crayon::num_ansi_colors' when loading 'VIBER'
#> Warning: replacing previous import 'cli::num_ansi_colors' by
#> 'crayon::num_ansi_colors' when loading 'easypar'
#> ✔ Loading ctree, 'Clone trees in cancer'. Support : <https://caravagn.github.io/ctree/>
#> Warning: replacing previous import 'crayon::%+%' by 'ggplot2::%+%' when loading
#> 'VIBER'
#> ✔ Loading VIBER, 'Variational inference for multivariate Binomial mixtures'. Support : <https://caravagn.github.io/VIBER/>
#> ✔ Loading lineaGT, 'Lineage inference from gene therapy'. Support : <https://caravagnalab.github.io/lineaGT/>
library(magrittr)
reticulate::use_condaenv("lineaGT")
```

The coverage dataset can be filtered calling the
[`filter_dataset()`](caravagnalab.github.io/lineaGT/reference/filter_dataset.md)
function.

``` r
data(cov.df.example)
data(vaf.df.example)
```

``` r
cov.example.filt = cov.df.example %>%
  filter_dataset(min_cov=5, min_frac=0.05)
#> ℹ Filtering the input dataset with minimum coverage 5 and minimum clusters frac…
#> ✔ Filtering the input dataset with minimum coverage 5 and minimum clusters frac…
#> 

cov.example.filt
#> # A tibble: 264 × 4
#>    IS    timepoints lineage coverage
#>    <chr> <chr>      <chr>      <int>
#>  1 IS100 t1         l1             0
#>  2 IS100 t2         l1           418
#>  3 IS100 t1         l2             0
#>  4 IS100 t2         l2            74
#>  5 IS101 t1         l1           502
#>  6 IS101 t2         l1           186
#>  7 IS101 t1         l2            62
#>  8 IS101 t2         l2           640
#>  9 IS11  t1         l1           128
#> 10 IS11  t2         l1           196
#> # ℹ 254 more rows
```

## Fitting the model

``` r
x = fit(
  cov.df = cov.example.filt,
  vaf.df = vaf.df.example,
  steps = 500,
  k_interval = c(5, 15),
  infer_growth = TRUE,
  timepoints_to_int = unlist(list("t1"=60, "t2"=150))
  )
#> ℹ Starting lineaGT model selection to retrieve the optimal number of clones
#> ✔ Starting lineaGT model selection to retrieve the optimal number of clones ...…
#> 
#> ℹ Fitting model to cluster ISs
#> ✔ Found 8 clones of ISs!
#> 
#> ℹ Fitting model to cluster mutations
#> ℹ Starting clustering of clone C0 mutations
#>  [ VIBER - variational fit ] 
#> 
#> ℹ Input n = 3, with k < 3. Dirichlet concentration α = 1e-06.
#> ℹ Starting clustering of clone C0 mutationsℹ Beta (a_0, b_0) = (1, 1); q_i = prior. Optimise: ε = 1e-10 or 5000 steps, r = 10 starts.
#> ℹ Starting clustering of clone C0 mutations
#> ✔ VIBER fit completed in 0.03 mins (status: converged)
#> ℹ Starting clustering of clone C0 mutations
#> ── [ VIBER ] My VIBER model n = 3 (w = 4 dimensions). Fit with k = 3 clusters. ─
#> ℹ Starting clustering of clone C0 mutations• Clusters: π = 67% [C1] and 33% [C3], with π > 0.
#> ℹ Starting clustering of clone C0 mutations• Binomials: θ = <0.09, 0.19, 0.01, 0> [C1] and <0.01, 0.36, 0.03, 0.4> [C3].
#> ℹ Starting clustering of clone C0 mutationsℹ Score(s): ELBO = -1461.853. Fit converged in 6 steps, ε = 1e-10.
#> ℹ Starting clustering of clone C0 mutations✔ Reduced to k = 2 (from 3) selecting VIBER cluster(s) with π > 0.166666666666667, and Binomial p > 0 in w > 0 dimension(s).
#> ℹ Starting clustering of clone C0 mutations✔ Starting clustering of clone C0 mutations ... done
#> 
#> ℹ Fitting model to cluster mutationsℹ Starting phylogeny inference of clone C0
#>  [ ctree ~ clone trees generator for C0 ] 
#> 
#> # A tibble: 3 × 8
#>   cluster   t1.l1 t2.l1  t1.l2 t2.l2 nMuts is.clonal is.driver
#>   <chr>     <dbl> <dbl>  <dbl> <dbl> <dbl> <lgl>     <lgl>    
#> 1 S1      0.0856  0.190 0      0         2 FALSE     TRUE     
#> 2 S2      0.00681 0.362 0.0288 0.396     1 FALSE     FALSE    
#> 3 C0      1       1     1      1         1 TRUE      FALSE
#> ✔ Trees per region 1, 2, 1, 1
#> ℹ Starting phylogeny inference of clone C0ℹ Total 2 tree structures - search is exahustive
#> ℹ Starting phylogeny inference of clone C0
#> ℹ Starting phylogeny inference of clone C0── Ranking trees 
#> ℹ Starting phylogeny inference of clone C0✔ 2  trees with non-zero score, storing 2
#> ℹ Starting phylogeny inference of clone C0✔ Starting phylogeny inference of clone C0 ... done
#> 
#> ℹ Fitting model to cluster mutationsℹ Starting clustering of clone C1 mutations
#>  [ VIBER - variational fit ] 
#> 
#> ℹ Input n = 8, with k < 8. Dirichlet concentration α = 1e-06.
#> ℹ Starting clustering of clone C1 mutationsℹ Beta (a_0, b_0) = (1, 1); q_i = prior. Optimise: ε = 1e-10 or 5000 steps, r = 10 starts.
#> ℹ Starting clustering of clone C1 mutations
#> ✔ VIBER fit completed in 0.04 mins (status: converged)
#> ℹ Starting clustering of clone C1 mutations
#> ── [ VIBER ] My VIBER model n = 8 (w = 4 dimensions). Fit with k = 8 clusters. ─
#> ℹ Starting clustering of clone C1 mutations• Clusters: π = 25% [C1], 25% [C4], 25% [C8], 13% [C3], and 13% [C5], with π >
#> 0.
#> ℹ Starting clustering of clone C1 mutations• Binomials: θ = <0.18, 0, 0.01, 0> [C1], <0, 0.08, 0.01, 0.15> [C4], <0.31, 0,
#> 0.01, 0.27> [C8], <0, 0.01, 0.02, 0> [C3], and <0, 0.28, 0.04, 0> [C5].
#> ℹ Starting clustering of clone C1 mutationsℹ Score(s): ELBO = -4566.431. Fit converged in 7 steps, ε = 1e-10.
#> ℹ Starting clustering of clone C1 mutations✔ Reduced to k = 5 (from 8) selecting VIBER cluster(s) with π > 0.0625, and Binomial p > 0 in w > 0 dimension(s).
#> ℹ Starting clustering of clone C1 mutations✔ Starting clustering of clone C1 mutations ... done
#> 
#> ℹ Fitting model to cluster mutationsℹ Starting phylogeny inference of clone C1
#>  [ ctree ~ clone trees generator for C1 ] 
#> 
#> # A tibble: 6 × 8
#>   cluster   t1.l1   t2.l1   t1.l2 t2.l2 nMuts is.clonal is.driver
#>   <chr>     <dbl>   <dbl>   <dbl> <dbl> <dbl> <lgl>     <lgl>    
#> 1 S1      0.184   0.00279 0       0         2 FALSE     FALSE    
#> 2 S2      0       0       0       0         1 FALSE     FALSE    
#> 3 S3      0.00162 0.0756  0.0122  0.147     2 FALSE     TRUE     
#> 4 S4      0.00405 0.283   0       0         1 FALSE     FALSE    
#> 5 S5      0.314   0.00281 0.00862 0.271     2 FALSE     FALSE    
#> 6 C1      1       1       1       1         1 TRUE      FALSE
#> ✔ Trees per region 2, 2, 1, 2
#> ℹ Starting phylogeny inference of clone C1ℹ Total 6 tree structures - search is exahustive
#> ℹ Starting phylogeny inference of clone C1✖ Starting phylogeny inference of clone C1 ... failed
#> 
#> ℹ Fitting model to cluster mutations
#> <subscriptOutOfBoundsError in model[var, ]: subscript out of bounds>
#> ℹ Starting clustering of clone C4 mutations
#>  [ VIBER - variational fit ] 
#> 
#> ℹ Input n = 6, with k < 6. Dirichlet concentration α = 1e-06.
#> ℹ Starting clustering of clone C4 mutationsℹ Beta (a_0, b_0) = (1, 1); q_i = prior. Optimise: ε = 1e-10 or 5000 steps, r = 10 starts.
#> ℹ Starting clustering of clone C4 mutations
#> ✔ VIBER fit completed in 0.04 mins (status: converged)
#> ℹ Starting clustering of clone C4 mutations
#> ── [ VIBER ] My VIBER model n = 6 (w = 4 dimensions). Fit with k = 6 clusters. ─
#> ℹ Starting clustering of clone C4 mutations• Clusters: π = 50% [C2], 17% [C1], 17% [C3], and 17% [C5], with π > 0.
#> ℹ Starting clustering of clone C4 mutations• Binomials: θ = <0.15, 0, 0.01, 0> [C2], <0, 0, 0.02, 0.29> [C1], <0.19, 0.22,
#> 0.3, 0.2> [C3], and <0.36, 0, 0.02, 0.28> [C5].
#> ℹ Starting clustering of clone C4 mutationsℹ Score(s): ELBO = -3594.675. Fit converged in 6 steps, ε = 1e-10.
#> ℹ Starting clustering of clone C4 mutations✔ Reduced to k = 4 (from 6) selecting VIBER cluster(s) with π > 0.0833333333333333, and Binomial p > 0 in w > 0 dimension(s).
#> ℹ Starting clustering of clone C4 mutations✔ Starting clustering of clone C4 mutations ... done
#> 
#> ℹ Fitting model to cluster mutationsℹ Starting phylogeny inference of clone C4
#>  [ ctree ~ clone trees generator for C4 ] 
#> 
#> # A tibble: 5 × 8
#>   cluster t1.l1   t2.l1  t1.l2 t2.l2 nMuts is.clonal is.driver
#>   <chr>   <dbl>   <dbl>  <dbl> <dbl> <dbl> <lgl>     <lgl>    
#> 1 S1      0     0       0.0186 0.292     1 FALSE     TRUE     
#> 2 S2      0.146 0.00159 0      0         3 FALSE     FALSE    
#> 3 S3      0.192 0.222   0.296  0.198     1 FALSE     FALSE    
#> 4 S4      0.358 0.00474 0.0170 0.275     1 FALSE     FALSE    
#> 5 C4      1     1       1      1         1 TRUE      FALSE
#> ✔ Trees per region 6, 1, 6, 5
#> ℹ Starting phylogeny inference of clone C4ℹ Total 54 tree structures - search is exahustive
#> ℹ Starting phylogeny inference of clone C4
#> ℹ Starting phylogeny inference of clone C4── Ranking trees 
#> ℹ Starting phylogeny inference of clone C4✔ 33  trees with non-zero score, storing 33
#> ℹ Starting phylogeny inference of clone C4✔ Starting phylogeny inference of clone C4 ... done
#> 
#> ℹ Fitting model to cluster mutationsℹ Starting clustering of clone C7 mutations
#>  [ VIBER - variational fit ] 
#> 
#> ℹ Input n = 2, with k < 2. Dirichlet concentration α = 1e-06.
#> ℹ Starting clustering of clone C7 mutationsℹ Beta (a_0, b_0) = (1, 1); q_i = prior. Optimise: ε = 1e-10 or 5000 steps, r = 10 starts.
#> ℹ Starting clustering of clone C7 mutations
#> ✔ VIBER fit completed in 0.03 mins (status: converged)
#> ℹ Starting clustering of clone C7 mutations
#> ── [ VIBER ] My VIBER model n = 2 (w = 4 dimensions). Fit with k = 2 clusters. ─
#> ℹ Starting clustering of clone C7 mutations• Clusters: π = 50% [C1] and 50% [C2], with π > 0.
#> ℹ Starting clustering of clone C7 mutations• Binomials: θ = <0.4, 0, 0.31, 0.35> [C1] and <0.01, 0.11, 0, 0> [C2].
#> ℹ Starting clustering of clone C7 mutationsℹ Score(s): ELBO = -1584.512. Fit converged in 5 steps, ε = 1e-10.
#> ℹ Starting clustering of clone C7 mutations✔ Starting clustering of clone C7 mutations ... done
#> 
#> ℹ Fitting model to cluster mutationsℹ Starting phylogeny inference of clone C7
#>  [ ctree ~ clone trees generator for C7 ] 
#> 
#> # A tibble: 3 × 8
#>   cluster  t1.l1   t2.l1 t1.l2 t2.l2 nMuts is.clonal is.driver
#>   <chr>    <dbl>   <dbl> <dbl> <dbl> <dbl> <lgl>     <lgl>    
#> 1 S1      0.396  0.00234 0.308 0.348     1 FALSE     FALSE    
#> 2 S2      0.0104 0.112   0     0         1 FALSE     TRUE     
#> 3 C7      1      1       1     1         1 TRUE      FALSE
#> ✔ Trees per region 2, 1, 1, 1
#> ℹ Starting phylogeny inference of clone C7ℹ Total 2 tree structures - search is exahustive
#> ℹ Starting phylogeny inference of clone C7
#> ℹ Starting phylogeny inference of clone C7── Ranking trees 
#> ℹ Starting phylogeny inference of clone C7✔ 2  trees with non-zero score, storing 2
#> ℹ Starting phylogeny inference of clone C7✔ Starting phylogeny inference of clone C7 ... done
#> 
#> ℹ Fitting model to cluster mutations✔ Fitting model to cluster mutations ... done
#> 
#> ℹ Fitting model to estimate population growth rates
#> ℹ Starting growth models inference of clone C0
#> [1] "C0"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  1.0000,   1.0000],
#>         [144.9737,  32.1573],
#>         [240.2606, 373.0132]], dtype=torch.float64)
#> [1] "C0.S1"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[ 0.0000,  0.0000],
#>         [12.4154,  0.0000],
#>         [45.6714,  0.0000]], dtype=torch.float64)
#> [1] "C0.S2"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  0.0000,   0.0000],
#>         [  0.9877,   0.9275],
#>         [ 87.0084, 147.8174]], dtype=torch.float64)
#> ✔ Starting growth models inference of clone C0 ... done
#> 
#> ℹ Fitting model to estimate population growth ratesℹ Starting growth models inference of clone C1
#> [1] "C1"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  1.0000,   1.0000],
#>         [245.1873,  22.7707],
#>         [177.3688, 289.0849]], dtype=torch.float64)
#> [1] "C1.S1"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[ 0.0000,  0.0000],
#>         [45.0851,  0.0000],
#>         [ 0.4955,  0.0000]], dtype=torch.float64)
#> [1] "C1.S2"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[0., 0.],
#>         [0., 0.],
#>         [0., 0.]], dtype=torch.float64)
#> [1] "C1.S3"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[ 0.0000,  0.0000],
#>         [ 0.3962,  0.2778],
#>         [13.4135, 42.4975]], dtype=torch.float64)
#> [1] "C1.S4"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[ 0.0000,  0.0000],
#>         [ 0.9934,  0.0000],
#>         [50.2584,  0.0000]], dtype=torch.float64)
#> [1] "C1.S5"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[ 0.0000,  0.0000],
#>         [76.9098,  0.1963],
#>         [ 0.4990, 78.4075]], dtype=torch.float64)
#> ✔ Starting growth models inference of clone C1 ... done
#> 
#> ℹ Fitting model to estimate population growth ratesℹ Starting growth models inference of clone C2
#> [1] "C2"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  1.0000,   1.0000],
#>         [  1.1279,   1.1034],
#>         [546.8077, 388.3853]], dtype=torch.float64)
#> ✔ Starting growth models inference of clone C2 ... done
#> 
#> ℹ Fitting model to estimate population growth ratesℹ Starting growth models inference of clone C3
#> [1] "C3"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  1.0000,   1.0000],
#>         [ 91.9191, 245.4643],
#>         [109.1249, 751.2116]], dtype=torch.float64)
#> [1] "C3.S1"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[ 0.0000,  0.0000],
#>         [17.2780,  0.0000],
#>         [ 0.0000,  0.0000]], dtype=torch.float64)
#> ✔ Starting growth models inference of clone C3 ... done
#> 
#> ℹ Fitting model to estimate population growth ratesℹ Starting growth models inference of clone C4
#> [1] "C4"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  1.0000,   1.0000],
#>         [284.7689,  51.3752],
#>         [209.2327, 491.7685]], dtype=torch.float64)
#> [1] "C4.S1"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  0.0000,   0.0000],
#>         [  0.0000,   0.9543],
#>         [  0.0000, 143.3544]], dtype=torch.float64)
#> [1] "C4.S3"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[ 0.0000,  0.0000],
#>         [54.6822, 15.2270],
#>         [46.3942, 97.5620]], dtype=torch.float64)
#> [1] "C4.S4"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  0.0000,   0.0000],
#>         [102.0333,   0.8729],
#>         [  0.9927, 135.3893]], dtype=torch.float64)
#> [1] "C4.S2"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[0.0000e+00, 0.0000e+00],
#>         [1.4882e+01, 0.0000e+00],
#>         [1.5798e-03, 0.0000e+00]], dtype=torch.float64)
#> ✔ Starting growth models inference of clone C4 ... done
#> 
#> ℹ Fitting model to estimate population growth ratesℹ Starting growth models inference of clone C5
#> [1] "C5"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[1.0000e+00, 1.0000e+00],
#>         [4.6677e-01, 9.0327e-01],
#>         [5.5075e+02, 8.2758e+02]], dtype=torch.float64)
#> [1] "C5.S1"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  0.0000,   0.0000],
#>         [  0.0000,   0.0000],
#>         [ 88.3387, 158.3940]], dtype=torch.float64)
#> ✔ Starting growth models inference of clone C5 ... done
#> 
#> ℹ Fitting model to estimate population growth ratesℹ Starting growth models inference of clone C6
#> [1] "C6"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  1.0000,   1.0000],
#>         [329.5520,  17.2995],
#>         [ 15.8495,  38.2649]], dtype=torch.float64)
#> ✔ Starting growth models inference of clone C6 ... done
#> 
#> ℹ Fitting model to estimate population growth ratesℹ Starting growth models inference of clone C7
#> [1] "C7"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[  1.0000,   1.0000],
#>         [  0.4702,   0.7794],
#>         [425.7162, 197.8717]], dtype=torch.float64)
#> [1] "C7.S1"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[ 0.0000,  0.0000],
#>         [ 0.1861,  0.2401],
#>         [ 0.9975, 68.8124]], dtype=torch.float64)
#> [1] "C7.S2"
#> tensor([[  0],
#>         [ 60],
#>         [150]], dtype=torch.int32)
#> tensor([[0.0000e+00, 0.0000e+00],
#>         [4.9108e-03, 0.0000e+00],
#>         [4.7749e+01, 0.0000e+00]], dtype=torch.float64)
#> ✔ Starting growth models inference of clone C7 ... done
#> 
#> ℹ Fitting model to estimate population growth rates✔ Fitting model to estimate population growth rates ... done
```

Printing the fitted object information regarding the data:

- lineages and timpoints present in the data,

- number of integration sites,

- number of inferred clones of ISs, estimated via model selection on the
  input range of number of clusters,

- for each clone, the number of assigned ISs and the mean coverage, per
  timepoint and lineage.

``` r
x
#> ── [ lineaGT ]  ──────── Python: /usr/share/miniconda/envs/lineaGT/bin/python ──
#> → Lineages: l1 and l2.
#> → Timepoints: t1 and t2.
#> → Number of Insertion Sites: 66.
#> 
#> ── Optimal IS model with k = 8.
#> 
#>     C4 (19 ISs) : l1 [285, 209]; l2 [ 51, 492] 
#>     C1 (15 ISs) : l1 [245, 177]; l2 [ 23, 289] 
#>      C0 (6 ISs) : l1 [145, 240]; l2 [ 32, 373] 
#>      C2 (6 ISs) : l1 [  1, 547]; l2 [  1, 388] 
#>      C3 (6 ISs) : l1 [ 92, 109]; l2 [245, 751] 
#>      C5 (6 ISs) : l1 [  0, 551]; l2 [  1, 828] 
#>      C6 (4 ISs) : l1 [330,  16]; l2 [ 17,  38] 
#>      C7 (4 ISs) : l1 [  0, 426]; l2 [  1, 198]
```
