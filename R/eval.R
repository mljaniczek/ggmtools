#' Confusion metrics for edge selection at a threshold
#'
#' Compares a score matrix (e.g. posterior inclusion probabilities) with a
#' true adjacency matrix. An edge is predicted when its score is
#' \code{>= threshold}. A 0/1 (or logical) estimated adjacency matrix can be
#' passed as \code{score_mat} with the default threshold.
#'
#' @param score_mat Numeric p x p matrix of edge scores, or a vector of edge
#'   scores (e.g. an upper triangle) of the same length as \code{truth_mat}.
#' @param truth_mat p x p matrix of true edges (non-zero = edge), or a vector
#'   matching \code{score_mat}.
#' @param threshold Numeric threshold. Default 0.5.
#' @param upper_only Logical. If \code{TRUE} (default), only the upper
#'   triangle is used (each undirected edge once); otherwise all off-diagonal
#'   entries. Ignored for vectors.
#'
#' @details Edges with a missing score or truth are dropped.
#'
#' @return A list with \code{TP}, \code{FP}, \code{TN}, \code{FN}, \code{TPR}
#'   (sensitivity, power), \code{FPR}, \code{precision}, \code{F1},
#'   \code{FDR} (\code{1 - precision}) and \code{MCC} (Matthews correlation
#'   coefficient). Rates whose denominator is zero are \code{NA}; in
#'   particular \code{FDR} is \code{NA} when no edge is predicted and
#'   \code{MCC} is \code{NA} when any margin of the confusion table is 0.
#' @export
#'
#' @examples
#' truth <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
#' score <- matrix(c(0, .9, .6, .9, 0, .2, .6, .2, 0), 3, 3)
#' confusion_at_threshold(score, truth, threshold = 0.5)
#'
#' # A 0/1 estimated adjacency matrix
#' confusion_at_threshold(threshold_to_adjacency(score, 0.5), truth)$MCC
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
  FDR <- ratio(FP, TP + FP)
  # as.numeric() avoids integer overflow in the products
  MCC <- ratio(as.numeric(TP) * TN - as.numeric(FP) * FN,
               sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  list(TP = TP, FP = FP, TN = TN, FN = FN, TPR = TPR, FPR = FPR,
       precision = precision, F1 = F1, FDR = FDR, MCC = MCC)
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
#'   likely an edge), or a vector of edge scores.
#' @param truth_mat p x p matrix of true edges (non-zero = edge), or a vector
#'   matching \code{score_mat}.
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

#' Precision-recall curve and AUPR for edge recovery
#'
#' Computes the precision-recall (PR) curve of a score matrix against a true
#' adjacency matrix, using the upper triangle. As in [roc_auc()], the curve
#' steps only at distinct score values (one point per block of tied
#' scores), so the result does not depend on the order of the edges.
#'
#' Three summaries of the curve are available:
#' \describe{
#'   \item{\code{"trapezoid"}}{Trapezoid rule over the curve points, starting
#'     from precision 1 at recall 0. This matches the tie-corrected
#'     \code{compute_aupr()} of the DiNA simulation. A straight line across
#'     a large tied block is not achievable in PR space, so this is
#'     optimistic for sparse estimators that tie many edges at 0.}
#'   \item{\code{"average_precision"}}{Step-wise average precision,
#'     \eqn{\sum_b (R_b - R_{b-1}) P_b} over tie blocks \eqn{b}, as in
#'     scikit-learn. A tied block contributes at the precision of the whole
#'     block.}
#'   \item{\code{"interpolated"}}{Exact area under the curve obtained by
#'     interpolating true and false positives linearly within each tied
#'     block (Davis & Goadrich 2006; the continuous version, as in the
#'     \pkg{PRROC} package). Precision is not linear in recall within a
#'     block, so this differs from the trapezoid.}
#' }
#' Without ties, \code{"average_precision"} and \code{"interpolated"} are
#' close; with large tied blocks \code{"trapezoid"} can be much larger than
#' both.
#'
#' Edges with a missing score or truth are dropped.
#'
#' @inheritParams roc_auc
#' @param method Summary of the curve: \code{"trapezoid"} (default),
#'   \code{"average_precision"} or \code{"interpolated"}. See Details.
#'
#' @return A list with \code{recall} and \code{precision} (the curve, one
#'   point per distinct score, preceded by the point (0, 1)),
#'   \code{thresholds} (the score at each point; \code{Inf} for the first)
#'   and \code{auc}. If the truth has no edges or no non-edges, \code{auc} is
#'   \code{NA}.
#'
#' @references Davis, J. & Goadrich, M. (2006). The relationship between
#'   precision-recall and ROC curves. *Proceedings of the 23rd International
#'   Conference on Machine Learning*, 233-240.
#' @seealso [roc_auc()]
#' @export
#'
#' @examples
#' truth <- matrix(0, 5, 5); truth[1, 2] <- truth[2, 3] <- 1
#' truth <- truth + t(truth)
#' score <- matrix(0, 5, 5); score[1, 2] <- 0.9; score[3, 4] <- 0.9
#' score <- score + t(score)
#' pr_auc(score, truth)$auc
#' pr_auc(score, truth, method = "average_precision")$auc
#' pr_auc(score, truth, method = "interpolated")$auc
pr_auc <- function(score_mat, truth_mat,
                   method = c("trapezoid", "average_precision", "interpolated")) {
  method <- match.arg(method)
  sv <- edge_vectors(score_mat, truth_mat, upper_only = TRUE)
  s <- sv$score; y <- sv$truth
  n_pos <- sum(y)
  if (n_pos == 0 || n_pos == length(y)) {
    return(list(recall = c(0, 1), precision = c(1, NA_real_),
                thresholds = c(Inf, -Inf), auc = NA_real_))
  }
  ord <- order(s, decreasing = TRUE)
  s <- s[ord]; y <- y[ord]
  last <- c(s[-1] != s[-length(s)], TRUE)          # last row of each tie block
  TP <- cumsum(y)[last]; FP <- cumsum(!y)[last]
  recall <- c(0, TP / n_pos)
  precision <- c(1, TP / (TP + FP))
  auc <- switch(method,
    trapezoid = sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2),
    average_precision = sum(diff(recall) * precision[-1]),
    interpolated = pr_interpolated_area(c(0, TP), c(0, FP)) / n_pos)
  list(recall = recall, precision = precision, thresholds = c(Inf, s[last]),
       auc = auc)
}

# Area under precision as a function of TP when TP and FP grow linearly
# within each block from (TP[b - 1], FP[b - 1]) to (TP[b], FP[b]); divide by
# the number of positives to get the area in recall units. Within a block,
# precision(x) = (a + x) / (c + r x) for x in [0, dT], with a = TP[b - 1],
# c = TP[b - 1] + FP[b - 1] and r = 1 + dF / dT.
pr_interpolated_area <- function(TP, FP) {
  area <- 0
  for (b in seq_along(TP)[-1]) {
    dT <- TP[b] - TP[b - 1]
    if (dT == 0) next
    a <- TP[b - 1]; c <- TP[b - 1] + FP[b - 1]
    r <- 1 + (FP[b] - FP[b - 1]) / dT
    area <- area + if (c == 0) dT / r else
      dT / r + (a - c / r) * log((c + r * dT) / c) / r
  }
  area
}

#' Overlap of the top-k edges of two score vectors
#'
#' Counts how many edges are among the \code{k} highest-scoring edges of
#' both \code{score1} and \code{score2}, e.g. two methods' differential
#' scores. Ties are handled without depending on the order of the edges:
#' when the k-th place falls inside a block of tied scores, each edge of
#' that block is in the top k with probability (places left) / (block
#' size), and the result is the expected overlap under random
#' tie-breaking, \eqn{\sum_e P(e \in A_k) P(e \in B_k)}. Without ties at the
#' boundary this is the ordinary overlap.
#'
#' Sorting with \code{order()} and taking the first \code{k} breaks ties by
#' position. For sparse estimators, which tie most edges at 0, that inflates
#' the overlap of two methods once \code{k} exceeds their number of non-zero
#' edges, because both take the same low-index edges first.
#'
#' @param score1,score2 Numeric vectors of edge scores of the same length, or
#'   p x p matrices (the upper triangle is used). Edges with a missing score
#'   in either are dropped.
#' @param k Integer vector of list sizes.
#'
#' @return A numeric vector of the same length as \code{k} with the
#'   (expected) overlap. Values of \code{k} larger than the number of edges
#'   are truncated.
#' @export
#'
#' @examples
#' a <- c(0.9, 0.8, 0, 0, 0, 0)
#' b <- c(0, 0.7, 0.6, 0, 0, 0)
#' topk_overlap(a, b, k = 1:6)
topk_overlap <- function(score1, score2, k) {
  if (!is.null(dim(score1))) score1 <- as.matrix(score1)[upper.tri(as.matrix(score1))]
  if (!is.null(dim(score2))) score2 <- as.matrix(score2)[upper.tri(as.matrix(score2))]
  if (length(score1) != length(score2))
    stop("score1 and score2 must have the same number of edges")
  ok <- !is.na(score1) & !is.na(score2)
  score1 <- score1[ok]; score2 <- score2[ok]
  vapply(k, function(kk) sum(topk_prob(score1, kk) * topk_prob(score2, kk)),
         numeric(1))
}

# P(edge is in the top k) under uniformly random tie-breaking
topk_prob <- function(s, k) {
  n <- length(s)
  n_above <- n - rank(s, ties.method = "max")      # strictly higher scores
  n_block <- rank(s, ties.method = "max") - rank(s, ties.method = "min") + 1
  pmin(1, pmax(0, (k - n_above) / n_block))
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
  if (is.null(dim(score_mat)) && is.null(dim(truth_mat))) {
    if (length(score_mat) != length(truth_mat))
      stop("score and truth vectors must have the same length")
    s <- score_mat
    y <- truth_mat != 0
    ok <- !is.na(s) & !is.na(y)
    return(list(score = s[ok], truth = y[ok]))
  }
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
