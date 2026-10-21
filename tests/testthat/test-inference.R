# DiNA's .bayes_fdr_call() / bayes_fdr_call() and multiGGMr's bayes_fdr()
# (before ggmtools) stopped at the first failing position, inside tie blocks
old_bayes_fdr_call <- function(pip_vec, fdr_target) {
  ord <- order(pip_vec, decreasing = TRUE, na.last = NA)
  running_fdr <- cumsum(1 - pip_vec[ord]) / seq_along(ord)
  k_star <- sum(cummax(running_fdr) <= fdr_target)
  called <- logical(length(pip_vec))
  if (k_star > 0) called[ord[seq_len(k_star)]] <- TRUE
  called
}

# DiNA's pooled_perm_pvalues() (tcga_demo/utils.R), absolute epsilon
old_pooled_pvalues <- function(obs_vec, null_pool) {
  asc_pool <- sort(null_pool[is.finite(null_pool)])
  n_lt <- findInterval(obs_vec - .Machine$double.eps, asc_pool)
  (1 + length(asc_pool) - n_lt) / (1 + length(asc_pool))
}

test_that("bayes_fdr selects the largest set with mean(1 - prob) <= alpha", {
  sel <- bayes_fdr(c(0.6, 0.99, 0.9, 0.97), alpha = 0.05)
  expect_equal(as.logical(sel), c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(attr(sel, "n_selected"), 3L)
  expect_equal(attr(sel, "expected_fdr"), mean(c(0.01, 0.03, 0.1)))
  expect_equal(attr(sel, "cutoff"), 0.9)
  expect_equal(attr(sel, "alpha"), 0.05)
  expect_equal(as.logical(sel), old_bayes_fdr_call(c(0.6, 0.99, 0.9, 0.97), 0.05))

  none <- bayes_fdr(c(0.5, 0.4), alpha = 0.05)
  expect_false(any(none))
  expect_identical(attr(none, "expected_fdr"), 0)
  expect_true(is.na(attr(none, "cutoff")))
  expect_error(bayes_fdr(c(0.5, 1.2)), "prob must lie")
  expect_error(bayes_fdr(0.5, alpha = 2), "alpha")
})

test_that("bayes_fdr keeps tied probabilities together (changed from DiNA/multiGGMr)", {
  prob <- c(1, 1, 0.9, 0.9, 0.9)
  # running mean of 1 - prob: 0, 0, 0.033, 0.05, 0.06
  old <- old_bayes_fdr_call(prob, 0.05)
  expect_equal(old, c(TRUE, TRUE, TRUE, TRUE, FALSE))   # 2 of 3 tied 0.9s
  new <- bayes_fdr(prob, alpha = 0.05)
  expect_equal(as.logical(new), c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(attr(new, "cutoff"), 1)
  # Old selection depends on input order, new does not
  expect_false(identical(old_bayes_fdr_call(rev(prob), 0.05), rev(old)) &&
               identical(old_bayes_fdr_call(prob[c(3, 1, 4, 2, 5)], 0.05),
                         old[c(3, 1, 4, 2, 5)]) &&
               identical(old_bayes_fdr_call(prob[c(5, 4, 3, 2, 1)], 0.05),
                         old[5:1]))
  for (perm in list(5:1, c(3, 1, 4, 2, 5)))
    expect_equal(as.logical(bayes_fdr(prob[perm], 0.05)), as.logical(new)[perm])
  # A whole tie block that fits is selected
  expect_true(all(bayes_fdr(c(1, 1, 0.97, 0.97, 0.97), 0.05)))
})

test_that("bayes_fdr degenerate and NA input", {
  expect_false(any(bayes_fdr(numeric(0))))
  expect_false(any(bayes_fdr(c(NA_real_, NA_real_))))
  sel <- bayes_fdr(c(0.99, NA, 0.98, 0.1))
  expect_equal(as.logical(sel), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(attr(sel, "expected_fdr"), 0.015)
  expect_true(all(bayes_fdr(rep(1, 4))))
  expect_false(any(bayes_fdr(rep(0, 4))))
})

test_that("bayes_fdr on a matrix uses the upper triangle and returns a symmetric matrix", {
  P <- matrix(0, 4, 4, dimnames = list(letters[1:4], letters[1:4]))
  P[upper.tri(P)] <- c(0.99, 0.2, 0.98, 0.1, 0.3, 0.97)
  P <- P + t(P)
  sel <- bayes_fdr(P, alpha = 0.05)
  expect_true(isSymmetric(unclass(sel) * 1))
  expect_false(any(diag(sel)))
  expect_equal(sum(sel[upper.tri(sel)]), 3)
  expect_equal(rownames(sel), letters[1:4])
  expect_error(bayes_fdr(matrix(0.5, 2, 3)), "square")
})

test_that("pooled_perm_pvalues matches spikeyglass and DiNA without ties", {
  set.seed(1)
  p <- 8
  null <- replicate(20, abs(matrix(rnorm(p * p), p)), simplify = FALSE)
  obs <- abs(matrix(rnorm(p * p), p))
  ut <- upper.tri(obs)
  pp <- pooled_perm_pvalues(obs, null)
  pool <- unlist(lapply(null, function(m) m[ut]))
  expect_equal(pp$pvalue, vapply(obs[ut], function(o)
    (1 + sum(pool >= o)) / (1 + length(pool)), numeric(1)))
  expect_equal(pp$pvalue, old_pooled_pvalues(obs[ut], pool))
  expect_equal(pp$p_adj, stats::p.adjust(pp$pvalue, "BH"))
  expect_length(pp$pvalue, p * (p - 1) / 2)
  # Vector input, null as a B x m matrix (spikeyglass's call)
  D <- do.call(rbind, lapply(null, function(m) m[ut]))
  expect_equal(pooled_perm_pvalues(obs[ut], D)$pvalue, pp$pvalue)
})

test_that("pooled_perm_pvalues counts ties as exceedances (changed from DiNA)", {
  null <- c(1, 2, 3, 3, 3, 4)
  obs <- c(3, 0.5)
  new <- pooled_perm_pvalues(obs, null)$pvalue
  expect_equal(new, c((1 + 4) / 7, 1))
  # DiNA's absolute epsilon vanishes for |obs| >= 2: ties at 3 not counted
  expect_equal(old_pooled_pvalues(obs, null), c((1 + 1) / 7, 1))
  # Near-ties from rounding count too
  expect_equal(pooled_perm_pvalues(0.1 + 0.2, 0.3)$pvalue, 1)
})

test_that("pooled_perm_pvalues handles NA and empty nulls", {
  pp <- pooled_perm_pvalues(c(1, NA, 3), list(c(0.5, NA, Inf), c(2, 2.5)))
  expect_equal(pp$pvalue, c(3 / 4, NA, 1 / 4))   # pool: 0.5, 2, 2.5
  expect_equal(pp$p_adj, stats::p.adjust(pp$pvalue, "BH"))
  expect_warning(e <- pooled_perm_pvalues(c(1, 2), c(NA, NaN)), "no finite")
  expect_true(all(is.na(e$pvalue)))
  expect_error(pooled_perm_pvalues(matrix(1, 2, 3), 1), "square")
})
