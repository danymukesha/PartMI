# Infer a direct association network using Part Mutual Information

Constructs a network where edges represent direct associations between
variables. For each pair (i, j), PMI is computed conditioning on all
other variables (full-order conditioning) or on a specified subset.
Statistical significance is assessed via permutation testing.

## Usage

``` r
pmi_network(
  data,
  method = c("binning", "knn"),
  bins = 5,
  k = 3,
  threshold = 0.05,
  n_permutations = 100,
  adjust = c("fdr", "bonferroni", "none"),
  mc_cores = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  Numeric matrix or data.frame (rows = samples, columns = variables).
  Variable/column names are used as node labels.

- method:

  PMI estimation method: `"binning"` or `"knn"`.

- bins:

  Number of bins for binning estimator (default 5).

- k:

  Number of neighbors for kNN estimator (default 3).

- threshold:

  Significance threshold for p-values (default 0.05).

- n_permutations:

  Number of permutations for significance testing (default 100). Use at
  least 1000 for published results.

- adjust:

  Method for multiple testing correction: `"none"`, `"bonferroni"`, or
  `"fdr"` (default `"fdr"`).

- mc_cores:

  Number of cores for parallel computation (default 1). Windows users
  should use 1.

- verbose:

  Print progress messages (default TRUE).

- ...:

  Additional arguments passed to
  [`pmi`](https://danymukesha.github.io/PartMI/reference/pmi.md).

## Value

A list with components:

- adjacency:

  Adjacency matrix (binary, 0/1) of significant edges.

- pmi_matrix:

  Matrix of PMI values for all pairs.

- pvalue_matrix:

  Matrix of permutation p-values.

- pvalue_adjusted:

  Matrix of adjusted p-values.

## References

Zhao et al. (2016) PNAS 113(18):5130-5135.

## Examples

``` r
set.seed(42)
n <- 50
g <- rnorm(n)
x <- g + rnorm(n, 0, 0.5)
y <- g + rnorm(n, 0, 0.5)
z <- x + rnorm(n, 0, 0.3)
w <- rnorm(n)
dat <- data.frame(X = x, Y = y, Z = z, W = w)
net <- pmi_network(dat, n_permutations = 20, threshold = 0.1)
#> Computing PMI for 4 variables...
#>   Pair 1/6 (17%)
#>   Pair 2/6 (33%)
#>   Pair 3/6 (50%)
#>   Pair 4/6 (67%)
#>   Pair 5/6 (83%)
#>   Pair 6/6 (100%)
#> Found 1 edges at threshold 0.100 (fdr corrected).
print(net$adjacency)
#>   X Y Z W
#> X 0 0 1 0
#> Y 0 0 0 0
#> Z 1 0 0 0
#> W 0 0 0 0
```
