sym <- function(m) { m[lower.tri(m)] <- t(m)[lower.tri(m)]; m }

# DiNA's compute_aupr() (04_calc_metrics.R, Oct 2026), for comparison
dina_compute_aupr <- function(truth_vec, score_vec) {
  valid <- !is.na(truth_vec) & !is.na(score_vec)
  truth_vec <- truth_vec[valid]; score_vec <- score_vec[valid]
  n_pos <- sum(truth_vec == 1)
  if (n_pos == 0 || n_pos == length(truth_vec)) return(NA_real_)
  ord <- order(score_vec, decreasing = TRUE)
  s <- score_vec[ord]; t <- truth_vec[ord]
  last <- c(s[-1] != s[-length(s)], TRUE)
  TP <- cumsum(t == 1)[last]; FP <- cumsum(t == 0)[last]
  precision <- TP / (TP + FP); recall <- TP / n_pos
  recall <- c(0, recall); precision <- c(1, precision)
  sum(diff(recall) * (precision[-length(precision)] + precision[-1]) / 2)
}

# Davis-Goadrich interpolated PR area by numerical integration: within each
# tie block TP and FP grow linearly; integrate precision over recall.
numeric_interpolated <- function(s, y) {
  o <- order(s, decreasing = TRUE); s <- s[o]; y <- y[o]
  last <- c(s[-1] != s[-length(s)], TRUE)
  TP <- c(0, cumsum(y)[last]); FP <- c(0, cumsum(!y)[last])
  area <- 0
  for (b in seq_along(TP)[-1]) {
    dT <- TP[b] - TP[b - 1]; dF <- FP[b] - FP[b - 1]
    if (dT == 0) next
    f <- function(x) (TP[b - 1] + x) / (TP[b - 1] + FP[b - 1] + x * (1 + dF / dT))
    area <- area + stats::integrate(f, 0, dT, rel.tol = 1e-10)$value
  }
  area / sum(y)
}

test_that("pr_auc trapezoid equals DiNA's compute_aupr and is order invariant", {
  set.seed(1)
  p <- 25
  truth <- sym(matrix(rbinom(p * p, 1, 0.2), p, p))
  score <- sym(matrix(ifelse(runif(p * p) < 0.7, 0, round(runif(p * p), 1)), p, p))
  ut <- upper.tri(score)
  r <- pr_auc(score, truth)
  expect_equal(r$auc, dina_compute_aupr(truth[ut], score[ut]))
  perm <- sample(p)
  for (m in c("trapezoid", "average_precision", "interpolated"))
    expect_equal(pr_auc(score, truth, m)$auc,
                 pr_auc(score[perm, perm], truth[perm, perm], m)$auc)
  expect_equal(r$recall[1], 0)
  expect_equal(r$precision[1], 1)
  expect_equal(tail(r$recall, 1), 1)
  expect_equal(length(r$thresholds), length(r$recall))
})

test_that("pr_auc methods on a small tied example", {
  # scores: 0.9 (pos), 0.5 tie block of (pos, neg, neg), 0.1 (pos)
  s <- c(0.9, 0.5, 0.5, 0.5, 0.1)
  y <- c(1, 1, 0, 0, 1)
  # points after (0, 1): (1/3, 1), (2/3, 1/2), (1, 3/5)
  expect_equal(pr_auc(s, y)$auc,
               (1/3) * 1 + (1/3) * (1 + 0.5) / 2 + (1/3) * (0.5 + 0.6) / 2)
  expect_equal(pr_auc(s, y, "average_precision")$auc,
               (1/3) * 1 + (1/3) * 0.5 + (1/3) * 0.6)
  expect_equal(pr_auc(s, y, "interpolated")$auc, numeric_interpolated(s, y))
  # Without ties every method reduces to the usual quantities
  s2 <- c(0.9, 0.8, 0.7, 0.6)
  y2 <- c(1, 0, 1, 0)
  ap <- (1 + 2/3) / 2
  expect_equal(pr_auc(s2, y2, "average_precision")$auc, ap)
  expect_equal(pr_auc(s2, y2, "trapezoid")$auc, 0.5 + 0.5 * (0.5 + 2/3) / 2)
})

test_that("interpolated PR area matches numerical integration with many ties", {
  set.seed(4)
  y <- rbinom(60, 1, 0.3)
  s <- ifelse(runif(60) < 0.8, 0, runif(60) + y)       # 80% tied at 0
  t <- pr_auc(s, y)$auc
  i <- pr_auc(s, y, "interpolated")$auc
  expect_equal(i, numeric_interpolated(s, y), tolerance = 1e-8)
})

test_that("trapezoid is optimistic for a large tied block at 0", {
  # Two top-ranked true edges, then 98 edges tied at 0 with 3 true edges
  s <- c(0.9, 0.8, rep(0, 98))
  y <- c(1, 1, rep(1, 3), rep(0, 95))
  expect_equal(pr_auc(s, y)$auc, 0.4 + 0.6 * (1 + 0.05) / 2)
  expect_equal(pr_auc(s, y, "average_precision")$auc, 0.4 + 0.6 * 0.05)
  i <- pr_auc(s, y, "interpolated")$auc
  expect_equal(i, numeric_interpolated(s, y))
  expect_lt(i, 0.6)
  expect_gt(i, 0.4 + 0.6 * 0.05)
})

test_that("pr_auc handles degenerate truth and NA", {
  expect_true(is.na(pr_auc(c(.1, .2, .3), c(0, 0, 0))$auc))
  expect_true(is.na(pr_auc(c(.1, .2, .3), c(1, 1, 1))$auc))
  expect_equal(pr_auc(c(.5, .5), c(0, 1))$auc, 0.75)   # one block: (0,1)-(1,1/2)
  expect_equal(pr_auc(c(.9, NA, .1, .5), c(1, 1, 0, NA))$auc,
               pr_auc(c(.9, .1), c(1, 0))$auc)
  expect_equal(pr_auc(c(.9, NA, .1, .5), c(1, 1, 0, NA))$auc, 1)
  expect_error(pr_auc(1:3, c(0, 1)), "same length")
})

test_that("roc_auc accepts vectors", {
  s <- c(0.9, 0.5, 0.5, 0.1); y <- c(1, 1, 0, 0)
  expect_equal(roc_auc(s, y)$auc, (1 + 0.5 + 1 + 1) / 4)
})

test_that("confusion_at_threshold FDR and MCC match DiNA's compute_fdr/compute_mcc", {
  dina_fdr <- function(TP, FP) if (TP + FP == 0) NA_real_ else FP / (TP + FP)
  dina_mcc <- function(TP, FP, TN, FN) {
    num <- as.numeric(TP) * TN - as.numeric(FP) * FN
    den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    if (den == 0) NA_real_ else num / den
  }
  set.seed(5)
  p <- 20
  truth <- sym(matrix(rbinom(p * p, 1, 0.2), p, p))
  score <- sym(matrix(round(runif(p * p), 1), p, p))
  for (thr in c(0.3, 0.5, 0.9)) {
    cm <- confusion_at_threshold(score, truth, threshold = thr)
    expect_equal(cm$FDR, dina_fdr(cm$TP, cm$FP))
    expect_equal(cm$MCC, dina_mcc(cm$TP, cm$FP, cm$TN, cm$FN))
    expect_equal(cm$FDR, 1 - cm$precision)
  }
  # Binary estimated adjacency with the default threshold
  est <- threshold_to_adjacency(score, 0.5)
  expect_equal(confusion_at_threshold(est, truth)[1:4],
               confusion_at_threshold(score, truth, threshold = 0.5 + 1e-9)[1:4])
})

test_that("confusion_at_threshold degenerate cases and NA", {
  truth <- sym(matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE))
  none <- confusion_at_threshold(matrix(0, 3, 3), truth)
  expect_true(is.na(none$FDR))                 # nothing predicted
  expect_true(is.na(none$MCC))
  expect_true(is.na(confusion_at_threshold(matrix(1, 3, 3), matrix(0, 3, 3))$MCC))
  expect_equal(confusion_at_threshold(c(1, 1, 0), c(1, 1, 0))$MCC, 1)
  expect_equal(confusion_at_threshold(c(0, 0, 1, 1), c(1, 1, 0, 0))$MCC, -1)
  cm <- confusion_at_threshold(c(1, NA, 0, 1), c(1, 0, NA, 0))
  expect_equal(c(cm$TP, cm$FP, cm$TN, cm$FN), c(1, 1, 0, 0))
})

test_that("topk_overlap is tie-aware and order invariant", {
  a <- c(0.9, 0.8, 0, 0, 0, 0)
  b <- c(0, 0.7, 0.6, 0, 0, 0)
  # k = 2: A = {1, 2}; B = {2, 3}; overlap 1
  expect_equal(topk_overlap(a, b, 2), 1)
  # k = 3: A = {1, 2} + one of 4 tied zeros (each 1/4); B = {2, 3, 1 of 4 zeros}
  # edge 1: 1 * 1/4; edge 2: 1; edge 3: 1/4 * 1; edges 4-6: 1/4 * 1/4 each
  expect_equal(topk_overlap(a, b, 3), 1/4 + 1 + 1/4 + 3/16)
  expect_equal(topk_overlap(a, b, 6), 6)
  expect_equal(topk_overlap(a, b, 10), 6)
  expect_equal(topk_overlap(a, b, 0), 0)

  # The order-based overlap (DiNA's plot_topk_overlap) depends on edge order
  naive <- function(x, y, k) length(intersect(order(x, decreasing = TRUE)[1:k],
                                               order(y, decreasing = TRUE)[1:k]))
  set.seed(6)
  x <- ifelse(runif(200) < 0.9, 0, runif(200))
  y <- ifelse(runif(200) < 0.9, 0, runif(200))
  perm <- sample(200)
  expect_equal(topk_overlap(x, y, c(10, 50, 100)),
               topk_overlap(x[perm], y[perm], c(10, 50, 100)))
  expect_false(naive(x, y, 50) == naive(x[perm], y[perm], 50) &&
               naive(x, y, 100) == naive(x[perm], y[perm], 100))
  # Expected overlap with all-tied scores is k^2 / n
  expect_equal(topk_overlap(rep(0, 10), rep(0, 10), 5), 2.5)
})

test_that("topk_overlap accepts matrices and drops NA", {
  m1 <- sym(matrix(c(0, .9, .1, 0, 0, .5, 0, 0, 0), 3, 3, byrow = TRUE))
  m2 <- sym(matrix(c(0, .8, .2, 0, 0, .1, 0, 0, 0), 3, 3, byrow = TRUE))
  expect_equal(topk_overlap(m1, m2, 1:3), c(1, 1, 3))
  expect_equal(topk_overlap(c(.9, NA, .1), c(.9, .5, .1), 1), 1)
  expect_error(topk_overlap(1:3, 1:4, 1), "same number")
})
