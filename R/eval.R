#' Confusion metrics for edge selection at a threshold
#'
#' Compares a score matrix (e.g. posterior inclusion probabilities) with a
#' true adjacency matrix. An edge is predicted when its score is
#' \code{>= threshold}.
#'
#' @param score_mat Numeric p x p matrix of edge scores.
#' @param truth_mat p x p matrix of true edges (non-zero = edge).
#' @param threshold Numeric threshold. Default 0.5.
#' @param upper_only Logical. If \code{TRUE} (default), only the upper
#'   triangle is used (each undirected edge once); otherwise all off-diagonal
#'   entries.
#'
#' @return A list with \code{TP}, \code{FP}, \code{TN}, \code{FN}, \code{TPR}
#'   (sensitivity), \code{FPR}, \code{precision} and \code{F1}. Rates whose
#'   denominator is zero are \code{NA}.
#' @export
#'
#' @examples
#' truth <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
#' score <- matrix(c(0, .9, .6, .9, 0, .2, .6, .2, 0), 3, 3)
#' confusion_at_threshold(score, truth, threshold = 0.5)
confusion_at_threshold <- function(score_mat, truth_mat, threshold = 0.5,
                                   upper_only = TRUE) {
  sv <- edge_vectors(score_mat, truth_mat, upper_only)
  pred <- sv$score >= threshold
  y <- sv$truth
  TP <- sum(pred & y); FP <- sum(pred & !y)
  TN <- sum(!pred & !y); FN <- sum(!pred & y)
  ratio <- function(a, b) if (b > 0) a / b else NA_real_
  TPR <- ratio(TP, TP + FN)
  FPR <- ratio(FP, FP + TN)
  precision <- ratio(TP, TP + FP)
  F1 <- if (is.na(precision) || is.na(TPR)) NA_real_ else
    if (precision + TPR > 0) 2 * precision * TPR / (precision + TPR) else 0
  list(TP = TP, FP = FP, TN = TN, FN = FN, TPR = TPR, FPR = FPR,
       precision = precision, F1 = F1)
}

#' ROC curve and AUC for edge recovery
#'
#' Computes the ROC curve of a score matrix against a true adjacency matrix,
#' using the upper triangle. The curve steps only at distinct score values,
#' so tied scores are handled correctly: a block of tied edges gives one
#' diagonal segment, i.e. each tied positive/negative pair counts 1/2. The
#' AUC therefore equals the Mann-Whitney statistic and does not depend on the
#' order of the edges.
#'
#' Ties matter for sparse estimators, which give many edges a score of
#' exactly 0 (and PIPs that underflow to exactly 0 or 1). An implementation
#' that accumulates one edge at a time in sorted order breaks ties by matrix
#' position and can be off by several hundredths of AUC.
#'
#' Edges with a missing score or truth are dropped.
#'
#' @param score_mat Numeric p x p matrix of edge scores (higher = more
#'   likely an edge).
#' @param truth_mat p x p matrix of true edges (non-zero = edge).
#'
#' @return A list with \code{FPR} and \code{TPR} (the curve, from (0, 0) to
#'   (1, 1)), \code{thresholds} (the score at each point; \code{Inf} for the
#'   origin) and \code{auc}. If the truth has no edges or no non-edges, the
#'   curve is the diagonal and \code{auc} is \code{NA}.
#' @export
#'
#' @examples
#' truth <- matrix(0, 5, 5); truth[1, 2] <- truth[2, 3] <- 1
#' truth <- truth + t(truth)
#' score <- matrix(0, 5, 5); score[1, 2] <- 0.9; score[3, 4] <- 0.9
#' score <- score + t(score)
#' roc_auc(score, truth)$auc   # ties at 0.9 and at 0 count 1/2
roc_auc <- function(score_mat, truth_mat) {
  sv <- edge_vectors(score_mat, truth_mat, upper_only = TRUE)
  s <- sv$score; y <- sv$truth
  n_pos <- sum(y); n_neg <- sum(!y)
  if (n_pos == 0 || n_neg == 0) {
    return(list(FPR = c(0, 1), TPR = c(0, 1), thresholds = c(Inf, -Inf),
                auc = NA_real_))
  }
  ord <- order(s, decreasing = TRUE)
  s <- s[ord]; y <- y[ord]
  last <- c(s[-1] != s[-length(s)], TRUE)          # last row of each tie block
  TPR <- c(0, cumsum(y)[last] / n_pos)
  FPR <- c(0, cumsum(!y)[last] / n_neg)
  auc <- sum(diff(FPR) * (TPR[-1] + TPR[-length(TPR)]) / 2)
  list(FPR = FPR, TPR = TPR, thresholds = c(Inf, s[last]), auc = auc)
}

#' Plot an ROC curve
#'
#' @param roc_obj A list from [roc_auc()]. (Older lists with an \code{AUC}
#'   element instead of \code{auc} are accepted.)
#' @param main Plot title.
#' @param add Logical. If \code{TRUE}, add the curve to the current plot.
#' @param ... Passed to \code{\link[graphics]{plot}} or
#'   \code{\link[graphics]{lines}}.
#'
#' @return \code{roc_obj}, invisibly.
#' @export
#'
#' @examples
#' truth <- matrix(rbinom(100, 1, 0.3), 10, 10); truth <- truth + t(truth)
#' plot_roc(roc_auc(matrix(runif(100), 10, 10), truth))
plot_roc <- function(roc_obj, main = "ROC curve", add = FALSE, ...) {
  auc <- roc_obj$auc %||% roc_obj$AUC
  if (add) {
    graphics::lines(roc_obj$FPR, roc_obj$TPR, ...)
  } else {
    graphics::plot(roc_obj$FPR, roc_obj$TPR, type = "l",
                   xlab = "False positive rate", ylab = "True positive rate",
                   main = main, xlim = c(0, 1), ylim = c(0, 1), ...)
    graphics::abline(a = 0, b = 1, lty = 2, col = "grey50")
    graphics::legend("bottomright", legend = sprintf("AUC = %.3f", auc),
                     bty = "n")
  }
  invisible(roc_obj)
}

edge_vectors <- function(score_mat, truth_mat, upper_only) {
  score_mat <- as.matrix(score_mat)
  truth_mat <- as.matrix(truth_mat)
  if (!identical(dim(score_mat), dim(truth_mat)))
    stop("score_mat and truth_mat must have the same dimensions")
  if (upper_only) {
    idx <- upper.tri(score_mat)
  } else {
    idx <- matrix(TRUE, nrow(score_mat), ncol(score_mat))
    diag(idx) <- FALSE
  }
  s <- score_mat[idx]
  y <- truth_mat[idx] != 0
  ok <- !is.na(s) & !is.na(y)
  list(score = s[ok], truth = y[ok])
}
