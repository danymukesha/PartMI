#' Estimate Mutual Information (MI) between two variables
#'
#' Computes the mutual information I(X; Y) using either a binning-based
#' histogram estimator or a k-nearest-neighbor (KSG) estimator.
#'
#' @param x Numeric vector (first variable).
#' @param y Numeric vector (second variable).
#' @param method Estimation method: \code{"binning"} (default) or \code{"knn"}.
#' @param bins Number of bins per dimension for the binning method (default 10).
#' @param k Number of nearest neighbors for the kNN method (default 3).
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric MI value in nats (natural log units).
#'
#' @references
#' Kraskov, Stögbauer & Grassberger (2004) Phys. Rev. E 69:066138.
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(100)
#' y <- x + rnorm(100, 0, 0.5)
#' mi(x, y)
#' mi(x, y, method = "knn", k = 5)
#'
#' @export
mi <- function(x, y, method = c("binning", "knn"), bins = 10, k = 3, ...) {
  method <- match.arg(method)
  switch(method,
    binning = .mi_binning(x, y, bins = bins),
    knn     = .mi_knn(x, y, k = k)
  )
}

#' Binning-based MI estimator
#' @keywords internal
.mi_binning <- function(x, y, bins = 10) {
  stopifnot(length(x) == length(y), length(x) > 0)
  # Remove NAs
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 2) return(0)

  xd <- .discretize(x, bins)
  yd <- .discretize(y, bins)

  tab <- table(xd, yd)
  p_xy <- tab / n
  p_x  <- rowSums(p_xy)
  p_y  <- colSums(p_xy)

  outer_p <- outer(p_x, p_y, "*")
  # Avoid division by zero
  mask <- p_xy > 0 & outer_p > 0
  if (!any(mask)) return(0)
  sum(p_xy[mask] * log(p_xy[mask] / outer_p[mask]))
}

#' kNN-based MI estimator (KSG algorithm 1)
#' @keywords internal
.mi_knn <- function(x, y, k = 3) {
  stopifnot(length(x) == length(y), length(x) > 0)
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n <= k + 1) {
    warning("Too few samples for kNN MI estimation")
    return(0)
  }

  xs <- .scale_unit(x)
  ys <- .scale_unit(y)

  # Precompute pairwise Chebyshev distances
  d_xy <- .chebyshev_dist(cbind(xs, ys))

  mi_val <- 0
  for (i in seq_len(n)) {
    d_i <- d_xy[i, -i]
    eps_i <- sort(d_i)[k]

    nx <- sum(abs(xs - xs[i]) < eps_i) - 1
    ny <- sum(abs(ys - ys[i]) < eps_i) - 1
    mi_val <- mi_val + .digamma(nx + 1) + .digamma(ny + 1)
  }

  .digamma(k) + .digamma(n) - mi_val / n
}
