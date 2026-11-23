#' Partial correlations from a joint graphical lasso fit
#'
#' [pcor_estimate()] method for fits from \code{JGL::JGL()} (class
#' \code{"jgl"}; Danaher, Wang & Witten 2014): converts each group's
#' estimated precision matrix with [precision_to_pcor()]. It only reads
#' \code{fit$theta}, so ggmtools does not depend on \pkg{JGL}. (JGL was
#' archived on CRAN in July 2025; it can be installed from the CRAN archive
#' or with \code{remotes::install_github("cran/JGL")}.)
#'
#' @param fit A \code{"jgl"} object from \code{JGL::JGL()}, fitted with
#'   \code{return.whole.theta = TRUE}. With the default
#'   \code{return.whole.theta = FALSE}, JGL drops variables that are not
#'   connected to any other, and the method stops with an error.
#' @param varnames Optional variable names. Defaults to the dimnames of the
#'   estimates (JGL keeps the column names of the data), then
#'   \code{V1, ..., Vp}.
#' @param groupnames Optional group names. JGL does not store them, so they
#'   default to \code{Group_1, ..., Group_K}; pass \code{names(data_list)} to
#'   use the names of the data list.
#' @param ... Ignored.
#'
#' @return A named list of K symmetric p x p partial-correlation matrices,
#'   as described in [pcor_estimate()]. These are penalized point estimates:
#'   entries the penalty sets to zero are exactly 0.
#' @references Danaher, P., Wang, P. & Witten, D.M. (2014). The joint
#'   graphical lasso for inverse covariance estimation across multiple
#'   classes. *Journal of the Royal Statistical Society: Series B*, 76(2),
#'   373-397.
#' @seealso [ggm_bic()] for tuning the penalties.
#' @export
#'
#' @examples
#' # A minimal object with the structure JGL::JGL() returns
#' theta <- list(diag(3), matrix(c(2, -0.5, 0, -0.5, 2, 0, 0, 0, 1), 3, 3))
#' fit <- structure(list(theta = theta, connected = rep(TRUE, 3)), class = "jgl")
#' pcor_estimate(fit, varnames = c("a", "b", "c"), groupnames = c("A", "B"))
pcor_estimate.jgl <- function(fit, varnames = NULL, groupnames = NULL, ...) {
  theta <- fit$theta
  if (!is.list(theta) || length(theta) < 1) stop("fit$theta must be a non-empty list")
  p <- nrow(theta[[1]])
  if (!is.null(fit$connected) && length(fit$connected) != p)
    stop("fit omits unconnected variables; refit JGL() with return.whole.theta = TRUE")
  varnames <- varnames %||% rownames(theta[[1]]) %||% paste0("V", seq_len(p))
  if (length(varnames) != p) stop("varnames must have length p")
  out <- lapply(theta, function(th) {
    P <- precision_to_pcor(th)
    dimnames(P) <- list(varnames, varnames)
    P
  })
  names(out) <- default_group_names(groupnames, length(out))
  out
}
