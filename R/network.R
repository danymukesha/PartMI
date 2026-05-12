#' Infer a direct association network using Part Mutual Information
#'
#' Constructs a network where edges represent direct associations between
#' variables. For each pair (i, j), PMI is computed conditioning on all other
#' variables (full-order conditioning) or on a specified subset. Statistical
#' significance is assessed via permutation testing.
#'
#' @param data Numeric matrix or data.frame (rows = samples, columns = variables).
#'   Variable/column names are used as node labels.
#' @param method PMI estimation method: \code{"binning"} or \code{"knn"}.
#' @param bins Number of bins for binning estimator (default 5).
#' @param k Number of neighbors for kNN estimator (default 3).
#' @param threshold Significance threshold for p-values (default 0.05).
#' @param n_permutations Number of permutations for significance testing
#'   (default 100). Use at least 1000 for published results.
#' @param adjust Method for multiple testing correction:
#'   \code{"none"}, \code{"bonferroni"}, or \code{"fdr"} (default \code{"fdr"}).
#' @param mc_cores Number of cores for parallel computation (default 1).
#'   Windows users should use 1.
#' @param verbose Print progress messages (default TRUE).
#' @param ... Additional arguments passed to \code{\link{pmi}}.
#'
#' @return A list with components:
#'   \item{adjacency}{Adjacency matrix (binary, 0/1) of significant edges.}
#'   \item{pmi_matrix}{Matrix of PMI values for all pairs.}
#'   \item{pvalue_matrix}{Matrix of permutation p-values.}
#'   \item{pvalue_adjusted}{Matrix of adjusted p-values.}
#'
#' @references
#' Zhao et al. (2016) PNAS 113(18):5130-5135.
#'
#' @examples
#' set.seed(42)
#' n <- 50
#' g <- rnorm(n)
#' x <- g + rnorm(n, 0, 0.5)
#' y <- g + rnorm(n, 0, 0.5)
#' z <- x + rnorm(n, 0, 0.3)
#' w <- rnorm(n)
#' dat <- data.frame(X = x, Y = y, Z = z, W = w)
#' net <- pmi_network(dat, n_permutations = 20, threshold = 0.1)
#' print(net$adjacency)
#'
#' @export
pmi_network <- function(data, method = c("binning", "knn"),
                        bins = 5, k = 3,
                        threshold = 0.05,
                        n_permutations = 100,
                        adjust = c("fdr", "bonferroni", "none"),
                        mc_cores = 1,
                        verbose = TRUE, ...) {
  method <- match.arg(method)
  adjust <- match.arg(adjust)

  data <- as.matrix(data)
  # Remove constant columns
  vars <- apply(data, 2, sd, na.rm = TRUE)
  keep <- vars > .Machine$double.eps^0.5
  if (sum(keep) < 2) stop("Need at least 2 variables with non-zero variance")
  data <- data[, keep, drop = FALSE]

  p <- ncol(data)
  var_names <- colnames(data)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(p))

  if (verbose) message("Computing PMI for ", p, " variables...")

  # Allocate matrices
  pmi_mat <- matrix(0, p, p)
  pval_mat <- matrix(1, p, p)

  n_pairs <- p * (p - 1) / 2
  pair_idx <- 0

  for (i in seq_len(p - 1)) {
    for (j in (i + 1):p) {
      pair_idx <- pair_idx + 1
      if (verbose && pair_idx %% max(1, floor(n_pairs / 10)) == 0) {
        message(sprintf("  Pair %d/%d (%.0f%%)", pair_idx, n_pairs,
                        100 * pair_idx / n_pairs))
      }

      # Z = all other variables
      z_idx <- setdiff(seq_len(p), c(i, j))
      z_data <- data[, z_idx, drop = FALSE]

      # Observed PMI
      args <- list(x = data[, i], y = data[, j], z = z_data,
                   method = method, bins = bins, k = k)
      pmi_obs <- do.call(pmi, args)

      # Permutation test
      n_greater <- 0
      n_perm <- if (n_permutations > 0) n_permutations else 0
      if (n_perm > 0) {
        for (perm in seq_len(n_perm)) {
          x_perm <- sample(data[, i])
          args_perm <- list(x = x_perm, y = data[, j], z = z_data,
                            method = method, bins = bins, k = k)
          pmi_perm <- do.call(pmi, args_perm)
          if (pmi_perm >= pmi_obs) n_greater <- n_greater + 1
        }
      }

      pval <- if (n_perm > 0) (n_greater + 1) / (n_perm + 1) else NA

      pmi_mat[i, j] <- pmi_obs
      pmi_mat[j, i] <- pmi_obs
      pval_mat[i, j] <- pval
      pval_mat[j, i] <- pval
    }
  }

  # Multiple testing correction
  p_vals <- pval_mat[upper.tri(pval_mat)]
  p_adj <- switch(adjust,
    none      = p_vals,
    bonferroni = p.adjust(p_vals, method = "bonferroni"),
    fdr       = p.adjust(p_vals, method = "fdr")
  )

  pval_adj_mat <- pval_mat
  pval_adj_mat[upper.tri(pval_adj_mat)] <- p_adj
  pval_adj_mat <- pmin(pval_adj_mat, t(pval_adj_mat), na.rm = TRUE)

  # Build adjacency matrix
  adj_mat <- matrix(0, p, p)
  adj_mat[pval_adj_mat < threshold & upper.tri(adj_mat)] <- 1
  adj_mat <- adj_mat + t(adj_mat)

  dimnames(adj_mat) <- list(var_names, var_names)
  dimnames(pmi_mat) <- list(var_names, var_names)
  dimnames(pval_mat) <- list(var_names, var_names)
  dimnames(pval_adj_mat) <- list(var_names, var_names)

  result <- list(
    adjacency       = adj_mat,
    pmi_matrix      = pmi_mat,
    pvalue_matrix   = pval_mat,
    pvalue_adjusted = pval_adj_mat,
    threshold       = threshold,
    method          = method,
    n_permutations  = n_permutations,
    adjust_method   = adjust
  )

  class(result) <- "pmi_network"

  if (verbose) {
    n_edges <- sum(adj_mat[upper.tri(adj_mat)])
    message(sprintf("Found %d edges at threshold %.3f (%s corrected).",
                    n_edges, threshold, adjust))
  }

  result
}

#' Print method for pmi_network objects
#' @param x A \code{pmi_network} object.
#' @param ... Additional arguments.
#' @export
print.pmi_network <- function(x, ...) {
  cat("PMI Network\n")
  cat("  Method:           ", x$method, "\n")
  cat("  Threshold:        ", x$threshold, "\n")
  cat("  Adjusted p-value: ", x$adjust_method, "\n")
  cat("  Permutations:     ", x$n_permutations, "\n")
  n_nodes <- nrow(x$adjacency)
  n_edges <- sum(x$adjacency[upper.tri(x$adjacency)])
  cat("  Nodes:            ", n_nodes, "\n")
  cat("  Edges:            ", n_edges, "\n")
  cat("  Density:          ",
      format(2 * n_edges / (n_nodes * (n_nodes - 1)), digits = 4), "\n")
  invisible(x)
}
