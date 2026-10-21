#' Bayesian false discovery rate selection
#'
#' Selects the largest set of hypotheses whose posterior expected false
#' discovery rate is at most \code{alpha}, using the direct posterior
#' probability approach (Newton et al. 2004). Hypotheses are ranked by their
#' posterior probability \code{prob}; the expected FDR of the top \eqn{k} is
#' the mean of \eqn{1 - prob} over them.
#'
#' @param prob Numeric vector of posterior probabilities (e.g. posterior
#'   inclusion probabilities), or a symmetric matrix of edge probabilities
#'   (only the upper triangle is used).
#' @param alpha Target posterior expected FDR. Default 0.05.
#'
#' @details The guarantee is only as good as the calibration of \code{prob}.
#'   For differential edges, pass a probability of a real difference such as
#'   \eqn{P(|D| > c)} or \eqn{P(g_k \ne g_m)}, not a sign probability.
#'
#'   Tied probabilities are selected or rejected together: the cutoff is
#'   only placed between distinct values of \code{prob}, so the selection
#'   does not depend on the order of the input. If including a whole block
#'   of ties would push the expected FDR above \code{alpha}, none of the
#'   block is selected.
#'
#'   Missing probabilities are never selected and do not count towards the
#'   expected FDR.
#'
#' @return A logical vector, or a symmetric logical matrix with \code{FALSE}
#'   diagonal, marking the selected hypotheses. Attributes:
#'   \code{expected_fdr} (posterior expected FDR of the selected set, 0 if
#'   none), \code{n_selected}, \code{cutoff} (smallest selected probability,
#'   or \code{NA}) and \code{alpha}.
#'
#' @references Newton, M.A., Noueiry, A., Sarkar, D. & Ahlquist, P. (2004).
#'   Detecting differential gene expression with a semiparametric
#'   hierarchical mixture method. *Biostatistics*, 5(2), 155-176.
#' @seealso [pooled_perm_pvalues()] for permutation-based selection.
#' @export
#'
#' @examples
#' bayes_fdr(c(0.99, 0.97, 0.9, 0.6, 0.2), alpha = 0.05)
#'
#' # Tied probabilities are kept together
#' bayes_fdr(c(1, 1, 0.9, 0.9, 0.9), alpha = 0.05)
bayes_fdr <- function(prob, alpha = 0.05) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha < 0 || alpha > 1)
    stop("alpha must be a single number in [0, 1].")
  is_mat <- is.matrix(prob)
  if (is_mat) {
    if (nrow(prob) != ncol(prob)) stop("prob matrix must be square.")
    v <- prob[upper.tri(prob)]
  } else {
    v <- as.numeric(prob)
  }
  if (any(v < 0 | v > 1, na.rm = TRUE)) stop("prob must lie in [0, 1].")

  ord <- order(v, decreasing = TRUE, na.last = NA)
  vs <- v[ord]
  running <- cumsum(1 - vs) / seq_along(vs)
  last <- c(vs[-1] != vs[-length(vs)], TRUE)      # last row of each tie block
  ok <- which(last & cummax(running) <= alpha)
  k <- if (length(ok)) max(ok) else 0L
  sel <- logical(length(v))
  if (k > 0) sel[ord[seq_len(k)]] <- TRUE

  if (is_mat) {
    m <- matrix(FALSE, nrow(prob), ncol(prob), dimnames = dimnames(prob))
    m[upper.tri(m)] <- sel
    sel <- m | t(m)
  }
  attr(sel, "expected_fdr") <- if (k > 0) running[k] else 0
  attr(sel, "n_selected") <- as.integer(k)
  attr(sel, "cutoff") <- if (k > 0) vs[k] else NA_real_
  attr(sel, "alpha") <- alpha
  sel
}

#' Empirical p-values against a pooled permutation null
#'
#' Computes one-sided empirical p-values for observed edge scores (larger =
#' more extreme, e.g. \eqn{|\hat\theta^{(1)}_{ij} - \hat\theta^{(2)}_{ij}|})
#' against a null distribution pooled over all edges and permutations, and
#' Benjamini-Hochberg adjusts them.
#'
#' Pooling gives a p-value resolution of about \eqn{1 / (1 + B m)} for
#' \eqn{B} permutations and \eqn{m} edges, instead of \eqn{1 / (1 + B)} for
#' per-edge nulls, which is fine enough to survive a BH correction over many
#' edges. It assumes the per-edge null distributions are approximately
#' exchangeable, which is reasonable for label-permutation nulls of
#' comparable (e.g. standardized) scores.
#'
#' @param observed Numeric vector of observed scores, or a square matrix
#'   (the upper triangle is used, in \code{upper.tri()} order as in
#'   [edge_index()]).
#' @param null Null scores: a numeric vector or matrix, or a list of them.
#'   All finite values are pooled. If \code{observed} is a matrix, list
#'   elements with the same dimensions contribute their upper triangle.
#' @param tol Relative tolerance: a null value counts as at least as large
#'   as \eqn{x} if it is \eqn{> x - tol \cdot \max(1, |x|)}, so that ties
#'   up to rounding error count as exceedances.
#'
#' @return A list with \code{pvalue}, \eqn{(1 + \#\{null \ge obs\}) / (1 +
#'   n_{null})}, and \code{p_adj} (BH-adjusted), one per observed edge.
#'   Missing observed scores give \code{NA}. If the pool has no finite
#'   values, all p-values are \code{NA} with a warning.
#' @seealso [bayes_fdr()]
#' @export
#'
#' @examples
#' set.seed(1)
#' p <- 6
#' null <- replicate(50, abs(matrix(rnorm(p * p), p)), simplify = FALSE)
#' obs <- abs(matrix(rnorm(p * p), p)); obs[1, 2] <- 4
#' pp <- pooled_perm_pvalues(obs, null)
#' cbind(edge_index(p), pvalue = pp$pvalue, p_adj = pp$p_adj)[1:3, ]
pooled_perm_pvalues <- function(observed, null, tol = sqrt(.Machine$double.eps)) {
  ut <- NULL
  if (is.matrix(observed)) {
    if (nrow(observed) != ncol(observed)) stop("observed matrix must be square")
    ut <- upper.tri(observed)
    obs <- as.numeric(observed[ut])
  } else {
    obs <- as.numeric(observed)
  }
  if (!is.list(null)) null <- list(null)
  pool <- unlist(lapply(null, function(m) {
    if (!is.null(ut) && is.matrix(m) && identical(dim(m), dim(ut))) m[ut]
    else as.vector(m)
  }), use.names = FALSE)
  pool <- sort(as.numeric(pool[is.finite(pool)]))
  if (length(pool) == 0) {
    warning("null has no finite values; returning NA")
    na <- rep(NA_real_, length(obs))
    return(list(pvalue = na, p_adj = na))
  }
  n_lt <- findInterval(obs - tol * pmax(1, abs(obs)), pool)
  pvalue <- (1 + length(pool) - n_lt) / (1 + length(pool))
  list(pvalue = pvalue, p_adj = stats::p.adjust(pvalue, method = "BH"))
}
