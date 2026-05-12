# Discretize continuous data into bins

Converts continuous expression values to discrete bins using
quantile-based or equal-width binning. Useful as a preprocessing step
for the binning PMI estimator.

## Usage

``` r
discretize_data(data, bins = 10, method = c("quantile", "width"))
```

## Arguments

- data:

  Numeric matrix or data.frame.

- bins:

  Number of bins (default 10).

- method:

  Binning method: `"quantile"` (equal frequency, default) or `"width"`
  (equal width).

## Value

An integer matrix of bin indices.

## Examples

``` r
set.seed(42)
dat <- matrix(rnorm(200), nrow = 20, ncol = 10)
d <- discretize_data(dat, bins = 5)
```
