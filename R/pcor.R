#' Convert a precision matrix to partial correlations
#'
#' Computes \eqn{\rho_{ij} = -\Omega_{ij} / \sqrt{\Omega_{ii}\Omega_{jj}}},
#' with unit diagonal.
#'
#' @param Omega A square precision (inverse covariance) matrix.
#'
#' @return A symmetric matrix of partial correlations with the dimnames of
#'   \code{Omega}. If any diagonal entry of \code{Omega} is not finite and
#'   positive, a matrix of \code{NA} is returned with a warning.
#' @export
#'
#' @examples
#' Omega <- matrix(c(2, -0.8, 0, -0.8, 2, -0.5, 0, -0.5, 1.5), 3, 3)
#' precision_to_pcor(Omega)
precision_to_pcor <- function(Omega) {
  Omega <- as.matrix(Omega)
  if (nrow(Omega) != ncol(Omega)) stop("Omega must be square")
  d <- diag(Omega)
  if (any(!is.finite(d)) || any(d <= 0)) {
    warning("Omega has a non-positive or non-finite diagonal; returning NA")
    out <- matrix(NA_real_, nrow(Omega), ncol(Omega), dimnames = dimnames(Omega))
    return(out)
  }
  P <- -Omega / sqrt(outer(d, d))
  P <- (P + t(P)) / 2
  diag(P) <- 1
  P
}
