sym <- function(m) { m[lower.tri(m)] <- t(m)[lower.tri(m)]; m }

mann_whitney_auc <- function(score, truth) {
  r <- rank(score); P <- sum(truth); N <- sum(!truth)
  (sum(r[truth]) - P * (P + 1) / 2) / (P * N)
}

test_that("roc_auc equals the Mann-Whitney AUC, also with ties", {
  set.seed(1)
  p <- 30
  truth <- sym(matrix(rbinom(p * p, 1, 0.2), p, p))
  for (score in list(
    sym(matrix(runif(p * p), p, p)),                       # no ties
    sym(matrix(round(runif(p * p), 1), p, p)),             # many ties
    sym(matrix(ifelse(runif(p * p) < 0.8, 0, runif(p * p)), p, p))  # sparse
  )) {
    ut <- upper.tri(score)
    r <- roc_auc(score, truth)
    expect_equal(r$auc, mann_whitney_auc(score[ut], truth[ut] != 0))
    expect_equal(c(r$FPR[1], r$TPR[1]), c(0, 0))
    expect_equal(c(tail(r$FPR, 1), tail(r$TPR, 1)), c(1, 1))
    expect_true(all(diff(r$FPR) >= 0) && all(diff(r$TPR) >= 0))
  }
})

test_that("roc_auc does not depend on the order of the variables", {
  set.seed(2)
  p <- 25
  truth <- sym(matrix(rbinom(p * p, 1, 0.3), p, p))
  score <- sym(matrix(ifelse(runif(p * p) < 0.7, 0, 1), p, p))  # two tie blocks
  perm <- sample(p)
  expect_equal(roc_auc(score, truth)$auc,
               roc_auc(score[perm, perm], truth[perm, perm])$auc)
})

test_that("roc_auc handles degenerate truth", {
  r <- roc_auc(matrix(runif(16), 4, 4), matrix(0, 4, 4))
  expect_true(is.na(r$auc))
})

test_that("confusion_at_threshold counts upper-triangle edges", {
  truth <- sym(matrix(c(0, 1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, 0,
                        0, 0, 0, 0), 4, 4, byrow = TRUE))
  score <- sym(matrix(c(0, .9, .6, 0,
                        0, 0, .2, .7,
                        0, 0, 0, 0,
                        0, 0, 0, 0), 4, 4, byrow = TRUE))
  cm <- confusion_at_threshold(score, truth, 0.5)
  expect_equal(cm[c("TP", "FP", "TN", "FN")], list(TP = 1L, FP = 2L, TN = 2L, FN = 1L))
  expect_equal(cm$TPR, 0.5)
  expect_equal(cm$FPR, 0.5)
  expect_equal(cm$precision, 1 / 3)
  expect_equal(cm$F1, 2 * (1 / 3) * 0.5 / (1 / 3 + 0.5))
  expect_true(is.na(confusion_at_threshold(score, 0 * truth)$TPR))
  expect_equal(confusion_at_threshold(score, truth, 0.5, upper_only = FALSE)$TP, 2L)
})

test_that("plot_roc draws and accepts old AUC element", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  r <- roc_auc(sym(matrix(runif(36), 6, 6)), sym(matrix(rbinom(36, 1, .4), 6, 6)))
  expect_invisible(plot_roc(r))
  expect_invisible(plot_roc(list(FPR = r$FPR, TPR = r$TPR, AUC = r$auc), add = TRUE))
})
