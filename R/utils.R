#' Scale vector to \[0, 1\]
#' @param x Numeric vector
#' @return Numeric vector scaled to \[0, 1\]
.scale_unit <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (r[1] == r[2]) {
        return(rep(0.5, length(x)))
    }
    (x - r[1]) / (r[2] - r[1])
}

#' Discretize a numeric vector into integer bin labels
#' @param x Numeric vector
#' @param bins Number of bins
#' @return Integer vector of bin indices (1..bins)
.discretize <- function(x, bins = 10) {
    if (length(unique(x)) <= bins) {
        return(match(x, sort(unique(x))))
    }
    breaks <- quantile(x, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE)
    breaks <- unique(breaks)
    as.integer(cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE))
}

#' Digamma function alias
#' @param x Numeric vector
.digamma <- function(x) {
    digamma(x)
}

#' Compute pairwise Chebyshev (L-infinity) distance matrix
#' @param mat Numeric matrix (rows = observations).
#' @return Symmetric distance matrix.
.chebyshev_dist <- function(mat) {
    n <- nrow(mat)
    d <- matrix(0, n, n)
    for (i in seq_len(n - 1)) {
        d[i, ] <- apply(mat, 1, function(row) max(abs(row - mat[i, ])))
        d[, i] <- d[i, ]
    }
    d
}
