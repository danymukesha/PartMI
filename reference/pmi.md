# Estimate Part Mutual Information (PMI) between X and Y given Z

Computes PMI(X; Y \| Z), quantifying the direct association between X
and Y after removing the influence of Z. Equivalent to conditional
mutual information (CMI) as applied to direct network inference.

## Usage

``` r
pmi(x, y, z, method = c("binning", "knn"), bins = 5, k = 3, ...)
```

## Arguments

- x:

  Numeric vector (first variable).

- y:

  Numeric vector (second variable).

- z:

  Numeric vector, matrix, or data.frame of conditioning variable(s).
  Each column is treated as a separate conditioning variable.

- method:

  Estimation method: `"binning"` (default) or `"knn"`.

- bins:

  Number of bins per dimension for binning (default 5). Decrease for
  higher-dimensional Z to avoid data sparsity.

- k:

  Number of nearest neighbors for the kNN method (default 3).

- ...:

  Additional arguments (currently unused).

## Value

A numeric PMI value in nats.

## Details

PMI is defined as: \$\$PMI(X;Y\|Z) = \sum\_{x,y,z} P(x,y,z) \log
\frac{P(x,y\|z)}{P(x\|z)P(y\|z)}\$\$

Two estimators are provided:

- `"binning"`: histogram-based estimator (fast, good for small Z)

- `"knn"`: k-nearest-neighbor estimator (Frenzel-Pompe, robust for
  higher-dimensional Z)

## References

Zhao et al. (2016) PNAS 113(18):5130-5135. Frenzel & Pompe (2007) Phys.
Rev. Lett. 99:204101.

## Examples

``` r
set.seed(42)
x <- rnorm(100)
z <- rnorm(100)
y <- x + z + rnorm(100, 0, 0.5)
# PMI between x and y conditioning on z should be small
pmi(x, y, z)
#>         1 
#> 0.6292994 
# MI (without conditioning) will typically be larger
mi(x, y)
#> [1] 0.7356992
```
