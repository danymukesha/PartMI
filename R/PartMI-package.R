#' PartMI: Part Mutual Information for Direct Association Network Inference
#'
#' Implements Part Mutual Information (PMI) for inferring direct associations
#' in biological networks. PMI extends standard Mutual Information (MI) by
#' conditioning on background variables, removing indirect dependencies and
#' reducing false-positive edges in network reconstruction.
#'
#' @section Core functions:
#' \itemize{
#'   \item \code{\link{mi}} - Compute mutual information between two variables
#'   \item \code{\link{pmi}} - Compute part (conditional) mutual information
#'   \item \code{\link{pmi_network}} - Infer direct association network
#'   \item \code{\link{normalize_data}} - Normalize expression data
#'   \item \code{\link{discretize_data}} - Discretize continuous data
#' }
#'
#' @references
#' Zhao et al. (2016) "Part mutual information for quantifying direct
#' associations in networks." PNAS 113(18): 5130-5135.
#'
#' Frenzel & Pompe (2007) "Partial Mutual Information for Coupling Analysis
#' of Multivariate Time Series." Physical Review Letters 99(20): 204101.
#'
#' Kraskov, Stögbauer & Grassberger (2004) "Estimating mutual information."
#' Physical Review E 69(6): 066138.
#'
#' @keywords internal
"_PACKAGE"
NULL
