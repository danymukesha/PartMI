# Normalize expression data

Applies quantile normalization or z-score normalization to expression
matrices. Designed for gene expression data (rows = samples, columns =
genes).

## Usage

``` r
normalize_data(data, method = c("zscore", "rank", "minmax", "quantile"), ...)
```

## Arguments

- data:

  Numeric matrix or data.frame (samples x variables).

- method:

  Normalization method:

  - `"zscore"` (default): Center and scale each column to mean 0, sd 1.

  - `"rank"`: Transform to normal scores via rank-based inverse normal.

  - `"minmax"`: Scale each column to \[0, 1\].

  - `"quantile"`: Quantile normalize across columns (samples must be
    aligned; i.e., rows are matched observations).

- ...:

  Additional arguments (currently unused).

## Value

Normalized matrix with same dimensions.

## Examples

``` r
set.seed(42)
dat <- matrix(rexp(200), nrow = 20, ncol = 10)
colnames(dat) <- paste0("Gene", 1:10)
norm_z <- normalize_data(dat, "zscore")
norm_r <- normalize_data(dat, "rank")
norm_q <- normalize_data(dat, "quantile")
```
