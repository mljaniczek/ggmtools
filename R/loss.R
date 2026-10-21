#' Relative Frobenius error of a matrix estimate
#'
#' Computes \eqn{\|\hat M - M\|_F / \|M\|_F} over all entries, including the
#' diagonal.
#'
#' @param est_mat Estimated matrix (e.g. a precision matrix).
#' @param truth_mat True matrix of the same dimensions.
#'
#' @return A number; \code{NA} if \code{truth_mat} is all zero or either
#'   matrix has missing entries.
#' @seealso [kl_divergence()]
#' @export
#'
#' @examples
#' Omega <- diag(3); Omega[1, 2] <- Omega[2, 1] <- 0.4
#' relative_frobenius(diag(3), Omega)
relative_frobenius <- function(est_mat, truth_mat) {
  est_mat <- as.matrix(est_mat); truth_mat <- as.matrix(truth_mat)
  if (!identical(dim(est_mat), dim(truth_mat)))
    stop("est_mat and truth_mat must have the same dimensions")
  truth_norm <- sqrt(sum(truth_mat^2))
  if (is.na(truth_norm) || truth_norm == 0) return(NA_real_)
  sqrt(sum((est_mat - truth_mat)^2)) / truth_norm
}

#' Kullback-Leibler divergence of an estimated Gaussian graphical model
#'
#' Computes the Kullback-Leibler divergence from the true zero-mean Gaussian
#' \eqn{N(0, \Sigma)} to the estimated one \eqn{N(0, \hat\Theta^{-1})},
#' \deqn{KL = \frac{1}{2}\left[\mathrm{tr}(\Sigma\hat\Theta) -
#'   \log|\Sigma\hat\Theta| - p\right].}
#' Twice this quantity is Stein's (entropy) loss of \eqn{\hat\Theta}.
#'
#' @param theta_hat Estimated p x p precision matrix.
#' @param sigma_true True p x p covariance matrix.
#'
#' @return A non-negative number, or \code{NA} if \eqn{\Sigma\hat\Theta}
#'   has a zero or negative determinant (e.g. \code{theta_hat} is not positive
#'   definite), the inputs have missing values, or the computation fails.
#' @seealso [relative_frobenius()]
#' @export
#'
#' @examples
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' kl_divergence(solve(Sigma), Sigma)   # 0
#' kl_divergence(diag(2), Sigma)
kl_divergence <- function(theta_hat, sigma_true) {
  theta_hat <- as.matrix(theta_hat); sigma_true <- as.matrix(sigma_true)
  if (!identical(dim(theta_hat), dim(sigma_true)))
    stop("theta_hat and sigma_true must have the same dimensions")
  if (anyNA(theta_hat) || anyNA(sigma_true)) return(NA_real_)
  tryCatch({
    prod_mat <- sigma_true %*% theta_hat
    log_det <- determinant(prod_mat, logarithm = TRUE)
    if (log_det$sign <= 0 || !is.finite(log_det$modulus)) return(NA_real_)
    (sum(diag(prod_mat)) - as.numeric(log_det$modulus) - nrow(theta_hat)) / 2
  }, error = function(e) NA_real_)
}

#' BIC for multiple Gaussian graphical models
#'
#' A BIC-type criterion for tuning joint estimators of K precision matrices
#' (e.g. the penalties of the joint graphical lasso), summed over groups:
#' \deqn{\sum_k \left[\mathrm{tr}(S_k\hat\Theta_k) - \log|\hat\Theta_k| +
#'   \frac{\log n_k}{n_k} E_k\right],}
#' where \eqn{S_k} is the sample covariance of group \eqn{k} (denominator
#' \eqn{n_k - 1}), \eqn{n_k} its sample size and \eqn{E_k} the number of
#' entries of \eqn{\hat\Theta_k} with \eqn{|\hat\theta_{ij}| \ge}
#' \code{threshold}. \eqn{E_k} counts every such entry of the matrix,
#' including the diagonal and both triangles. This is the per-observation
#' form used in the Koyejo lab joint graphical lasso code
#' (\code{vBIC()}); smaller is better.
#'
#' @param data_list List of K data matrices (\eqn{n_k} x p).
#' @param theta_list List of K estimated p x p precision matrices.
#' @param threshold Entries with absolute value \code{>= threshold} count as
#'   non-zero. Default 0.001.
#'
#' @return A number; \code{NA} if any \eqn{\hat\Theta_k} has missing values
#'   or a zero or negative determinant, or the data have missing values.
#' @references Danaher, P., Wang, P. & Witten, D.M. (2014). The joint
#'   graphical lasso for inverse covariance estimation across multiple
#'   classes. *Journal of the Royal Statistical Society: Series B*, 76(2),
#'   373-397.
#' @export
#'
#' @examples
#' set.seed(1)
#' X <- list(matrix(rnorm(100), 25, 4), matrix(rnorm(80), 20, 4))
#' theta <- lapply(X, function(x) solve(cov(x)))
#' ggm_bic(X, theta)
ggm_bic <- function(data_list, theta_list, threshold = 0.001) {
  if (!is.list(data_list) || !is.list(theta_list) ||
      length(data_list) != length(theta_list) || length(data_list) < 1)
    stop("data_list and theta_list must be lists of the same, non-zero length")
  total <- 0
  for (k in seq_along(data_list)) {
    X <- as.matrix(data_list[[k]])
    theta <- as.matrix(theta_list[[k]])
    n <- nrow(X)
    if (anyNA(theta)) return(NA_real_)
    log_det <- determinant(theta, logarithm = TRUE)
    if (log_det$sign <= 0 || !is.finite(log_det$modulus)) return(NA_real_)
    E <- sum(abs(theta) >= threshold)
    total <- total + sum(diag(stats::cov(X) %*% theta)) -
      as.numeric(log_det$modulus) + log(n) * E / n
  }
  total
}
