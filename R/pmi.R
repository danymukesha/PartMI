#' Estimate Part Mutual Information (PMI) between X and Y given Z
#'
#' Computes PMI(X; Y | Z), quantifying the direct association between X and Y
#' after removing the influence of Z. Equivalent to conditional mutual information
#' (CMI) as applied to direct network inference.
#'
#' PMI is defined as:
#' \deqn{PMI(X;Y|Z) = \sum_{x,y,z} P(x,y,z) \log \frac{P(x,y|z)}{P(x|z)P(y|z)}}
#'
#' Two estimators are provided:
#' \itemize{
#'   \item \code{"binning"}: histogram-based estimator (fast, good for small Z)
#'   \item \code{"knn"}: k-nearest-neighbor estimator (Frenzel-Pompe, robust for
#'     higher-dimensional Z)
#' }
#'
#' @param x Numeric vector (first variable).
#' @param y Numeric vector (second variable).
#' @param z Numeric vector, matrix, or data.frame of conditioning variable(s).
#'   Each column is treated as a separate conditioning variable.
#' @param method Estimation method: \code{"binning"} (default) or \code{"knn"}.
#' @param bins Number of bins per dimension for binning (default 5). Decrease
#'   for higher-dimensional Z to avoid data sparsity.
#' @param k Number of nearest neighbors for the kNN method (default 3).
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric PMI value in nats.
#'
#' @references
#' Zhao et al. (2016) PNAS 113(18):5130-5135.
#' Frenzel & Pompe (2007) Phys. Rev. Lett. 99:204101.
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(100)
#' z <- rnorm(100)
#' y <- x + z + rnorm(100, 0, 0.5)
#' # PMI between x and y conditioning on z should be small
#' pmi(x, y, z)
#' # MI (without conditioning) will typically be larger
#' mi(x, y)
#'
#' @export
pmi <- function(x, y, z, method = c("binning", "knn"),
                bins = 5, k = 3, ...) {
  method <- match.arg(method)
  # Handle empty Z: PMI with no conditioning = MI
  if (is.null(z) || (is.vector(z) && length(z) == 0) ||
      (is.matrix(z) && ncol(z) == 0)) {
    return(mi(x, y, method = method, bins = bins, k = k))
  }
  switch(method,
    binning = .pmi_binning(x, y, z, bins = bins),
    knn     = .pmi_knn(x, y, z, k = k)
  )
}

#' PMI via multidimensional histogram
#' @keywords internal
.pmi_binning <- function(x, y, z, bins = 5) {
  stopifnot(length(x) == length(y))
  # Coerce z to matrix
  z <- as.matrix(z)
  stopifnot(nrow(z) == length(x))

  ok <- complete.cases(x, y, z)
  x <- x[ok]; y <- y[ok]; z <- z[ok, , drop = FALSE]
  n <- length(x)
  if (n < 2) return(0)

  # Determine bins per dimension for Z
  nz <- ncol(z)
  # Adjust bins to avoid empty cells: reduce bins for high-dim Z
  bins_z <- max(2, floor(bins / nz))

  xd <- .discretize(x, bins)
  yd <- .discretize(y, bins)
  zd <- apply(z, 2, .discretize, bins = bins_z)
  if (nz == 1) zd <- matrix(zd, ncol = 1)

  # Create a single factor from the Z discretization
  if (nz == 1) {
    z_factor <- zd[, 1]
  } else {
    # Combine multi-dim Z into one factor
    z_factor <- apply(zd, 1, paste, collapse = ":")
  }
  z_factor <- as.factor(z_factor)
  z_levels <- nlevels(z_factor)

  # Build 3D contingency table via looping over Z levels
  pmi_val <- 0
  for (zv in seq_len(z_levels)) {
    idx <- which(as.integer(z_factor) == zv)
    if (length(idx) < 2) next
    n_z <- length(idx)

    x_sub <- xd[idx]
    y_sub <- yd[idx]

    tab_xy <- table(factor(x_sub, levels = seq_len(bins)),
                    factor(y_sub, levels = seq_len(bins)))
    tab_x  <- rowSums(tab_xy) / n_z
    tab_y  <- colSums(tab_xy) / n_z
    p_xy   <- tab_xy / n_z

    p_z <- n_z / n

    for (i in seq_len(bins)) {
      for (j in seq_len(bins)) {
        pxy <- p_xy[i, j]
        px  <- tab_x[i]
        py  <- tab_y[j]
        if (pxy > 0 && px > 0 && py > 0) {
          pmi_val <- pmi_val + p_z * pxy * log(pxy / (px * py))
        }
      }
    }
  }

  pmi_val
}

#' PMI via kNN (Frenzel-Pompe estimator)
#'
#' Uses precomputed distance matrices for efficiency.
#' @keywords internal
.pmi_knn <- function(x, y, z, k = 3) {
  stopifnot(length(x) == length(y))
  z <- as.matrix(z)
  stopifnot(nrow(z) == length(x))

  ok <- complete.cases(x, y, z)
  x <- x[ok]; y <- y[ok]; z <- z[ok, , drop = FALSE]
  n <- length(x)
  if (n <= k + 1) {
    warning("Too few samples for kNN PMI estimation")
    return(0)
  }

  xs <- .scale_unit(x)
  ys <- .scale_unit(y)
  zs <- apply(z, 2, .scale_unit)

  # Precompute pairwise Chebyshev distance matrices
  # d_ab[i,j] = max(|a_i - a_j|, |b_i - b_j|)
  d_z  <- .chebyshev_dist(zs)
  d_xz <- .chebyshev_dist(cbind(xs, zs))
  d_yz <- .chebyshev_dist(cbind(ys, zs))
  d_full <- .chebyshev_dist(cbind(xs, ys, zs))

  pmi_val <- 0

  for (i in seq_len(n)) {
    # Distance to k-th NN in full space (exclude self)
    d_i <- d_full[i, -i]
    eps_i <- sort(d_i)[k]

    # Count neighbors within eps (exclude self)
    n_z  <- sum(d_z[i, -i]  < eps_i)
    n_xz <- sum(d_xz[i, -i] < eps_i)
    n_yz <- sum(d_yz[i, -i] < eps_i)

    pmi_val <- pmi_val + .digamma(n_z + 1) - .digamma(n_xz + 1) - .digamma(n_yz + 1)
  }

  .digamma(k) + pmi_val / n
}
