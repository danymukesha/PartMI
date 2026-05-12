# Estimate Mutual Information (MI) between two variables

Computes the mutual information I(X; Y) using either a binning-based
histogram estimator or a k-nearest-neighbor (KSG) estimator.

## Usage

``` r
mi(x, y, method = c("binning", "knn"), bins = 10, k = 3, ...)
```

## Arguments

- x:

  Numeric vector (first variable).

- y:

  Numeric vector (second variable).

- method:

  Estimation method: `"binning"` (default) or `"knn"`.

- bins:

  Number of bins per dimension for the binning method (default 10).

- k:

  Number of nearest neighbors for the kNN method (default 3).

- ...:

  Additional arguments (currently unused).

## Value

A numeric MI value in nats (natural log units).

## References

Kraskov, Stögbauer & Grassberger (2004) Phys. Rev. E 69:066138.

## Examples

``` r
set.seed(42)
x <- rnorm(100)
y <- x + rnorm(100, 0, 0.5)
mi(x, y)
#> [1] 1.011395
mi(x, y, method = "knn", k = 5)
#> [1] 0.9244922
```
