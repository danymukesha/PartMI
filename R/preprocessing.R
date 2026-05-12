#' Normalize expression data
#'
#' Applies quantile normalization or z-score normalization to expression
#' matrices. Designed for gene expression data (rows = samples, columns = genes).
#'
#' @param data Numeric matrix or data.frame (samples x variables).
#' @param method Normalization method:
#'   \itemize{
#'     \item \code{"zscore"} (default): Center and scale each column to mean 0, sd 1.
#'     \item \code{"rank"}: Transform to normal scores via rank-based inverse normal.
#'     \item \code{"minmax"}: Scale each column to \[0, 1\].
#'     \item \code{"quantile"}: Quantile normalize across columns (samples must be
#'       aligned; i.e., rows are matched observations).
#'   }
#' @param ... Additional arguments (currently unused).
#'
#' @return Normalized matrix with same dimensions.
#'
#' @examples
#' set.seed(42)
#' dat <- matrix(rexp(200), nrow = 20, ncol = 10)
#' colnames(dat) <- paste0("Gene", 1:10)
#' norm_z <- normalize_data(dat, "zscore")
#' norm_r <- normalize_data(dat, "rank")
#' norm_q <- normalize_data(dat, "quantile")
#'
#' @export
normalize_data <- function(data, method = c("zscore", "rank", "minmax", "quantile"),
                           ...) {
    method <- match.arg(method)
    data <- as.matrix(data)

    switch(method,
        zscore = {
            means <- colMeans(data, na.rm = TRUE)
            sds <- apply(data, 2, sd, na.rm = TRUE)
            sds[sds == 0] <- 1
            scale(data, center = means, scale = sds)
        },
        rank = {
            .rank_intnorm(data)
        },
        minmax = {
            apply(data, 2, function(x) {
                r <- range(x, na.rm = TRUE)
                if (r[1] == r[2]) {
                    return(rep(0.5, length(x)))
                }
                (x - r[1]) / (r[2] - r[1])
            })
        },
        quantile = {
            .quantile_norm(data)
        }
    )
}

#' Rank-based inverse normal transformation
#' @keywords internal
.rank_intnorm <- function(x) {
    apply(x, 2, function(col) {
        n <- sum(!is.na(col))
        if (n < 2) {
            return(col)
        }
        ranks <- rank(col, na.last = "keep", ties.method = "average")
        qnorm((ranks - 0.5) / n)
    })
}

#' Simple quantile normalization
#' @keywords internal
.quantile_norm <- function(x) {
    x_rank <- apply(x, 2, rank, ties.method = "average")
    x_sorted <- apply(x, 2, sort)
    q_mean <- rowMeans(x_sorted)
    # Map back
    apply(x_rank, 2, function(r) q_mean[round(r)])
}

#' Discretize continuous data into bins
#'
#' Converts continuous expression values to discrete bins using quantile-based
#' or equal-width binning. Useful as a preprocessing step for the binning
#' PMI estimator.
#'
#' @param data Numeric matrix or data.frame.
#' @param bins Number of bins (default 10).
#' @param method Binning method: \code{"quantile"} (equal frequency, default)
#'   or \code{"width"} (equal width).
#'
#' @return An integer matrix of bin indices.
#'
#' @examples
#' set.seed(42)
#' dat <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' d <- discretize_data(dat, bins = 5)
#'
#' @export
discretize_data <- function(data, bins = 10, method = c("quantile", "width")) {
    method <- match.arg(method)
    data <- as.matrix(data)

    result <- matrix(0, nrow(data), ncol(data))
    for (j in seq_len(ncol(data))) {
        result[, j] <- .discretize_col(data[, j], bins, method)
    }
    colnames(result) <- colnames(data)
    rownames(result) <- rownames(data)
    result
}

.discretize_col <- function(x, bins, method) {
    if (length(unique(x)) <= bins) {
        return(match(x, sort(unique(x))))
    }
    if (method == "quantile") {
        breaks <- quantile(x, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE)
    } else {
        r <- range(x, na.rm = TRUE)
        breaks <- seq(r[1], r[2], length.out = bins + 1)
    }
    breaks <- unique(breaks)
    as.integer(cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE))
}
