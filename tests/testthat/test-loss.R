# DiNA's compute_kl_divergence() (04_calc_metrics.R) and vBIC()
# (00_helper_functions.R), for comparison
dina_kl <- function(theta_hat, sigma_true) {
  prod_mat <- sigma_true %*% theta_hat
  log_det <- determinant(prod_mat, logarithm = TRUE)
  if (log_det$sign <= 0) return(NA_real_)
  sum(diag(prod_mat)) - log_det$modulus[1] - nrow(theta_hat)
}
dina_vBIC <- function(X, est_graph, thr = 0.001) {
  BIC_acc <- 0
  for (i in seq_along(X)) {
    data_num <- dim(X[[i]])[1]
    sample_cov <- cov(X[[i]], X[[i]])
    tr_sum <- sum(diag(sample_cov %*% est_graph$theta[[i]]))
    log_det <- determinant(est_graph$theta[[i]], logarithm = TRUE)$modulus[1][1]
    E <- sum(sum(abs(est_graph$theta[[i]]) >= thr))
    BIC_acc <- BIC_acc + (tr_sum - log_det) + (log(data_num) * E / data_num)
  }
  BIC_acc
}

test_that("relative_frobenius", {
  truth <- matrix(c(2, 1, 1, 2), 2, 2)
  expect_equal(relative_frobenius(truth, truth), 0)
  expect_equal(relative_frobenius(2 * truth, truth), 1)
  expect_equal(relative_frobenius(diag(2), truth), sqrt(1 + 1 + 1 + 1) / sqrt(10))
  expect_true(is.na(relative_frobenius(diag(2), matrix(0, 2, 2))))
  expect_true(is.na(relative_frobenius(matrix(c(1, NA, NA, 1), 2), truth)))
  expect_true(is.na(relative_frobenius(diag(2), matrix(c(1, NA, NA, 1), 2))))
  expect_error(relative_frobenius(diag(2), diag(3)), "same dimensions")
})

test_that("kl_divergence is the Gaussian KL, half of DiNA's value", {
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  expect_equal(kl_divergence(solve(Sigma), Sigma), 0)
  # Closed form KL(N(0, Sigma) || N(0, I)) = (tr(Sigma) - log|Sigma| - p) / 2
  expect_equal(kl_divergence(diag(2), Sigma), (2 - log(0.75) - 2) / 2)
  # DiNA's compute_kl_divergence() returned twice this (Stein's loss)
  set.seed(2)
  A <- matrix(rnorm(16), 4); S <- crossprod(A) + diag(4)
  Th <- solve(S + diag(0.3, 4))
  expect_equal(dina_kl(Th, S), 2 * kl_divergence(Th, S))
  expect_equal(dina_kl(diag(2), Sigma), 2 * kl_divergence(diag(2), Sigma))
  expect_gt(kl_divergence(Th, S), 0)
})

test_that("kl_divergence returns NA for invalid input", {
  Sigma <- diag(2)
  expect_true(is.na(kl_divergence(diag(c(1, -1)), Sigma)))      # negative det
  expect_true(is.na(kl_divergence(matrix(0, 2, 2), Sigma)))     # singular
  expect_true(is.na(kl_divergence(matrix(c(1, NA, NA, 1), 2), Sigma)))
  expect_error(kl_divergence(diag(2), diag(3)), "same dimensions")
})

test_that("ggm_bic equals DiNA's vBIC", {
  set.seed(1)
  X <- list(matrix(rnorm(100), 25, 4), matrix(rnorm(80), 20, 4))
  theta <- lapply(X, function(x) { th <- solve(cov(x)); th[abs(th) < 0.1] <- 0; th })
  expect_equal(ggm_bic(X, theta), dina_vBIC(X, list(theta = theta)))
  expect_equal(ggm_bic(X, theta, threshold = 0.2),
               dina_vBIC(X, list(theta = theta), thr = 0.2))
  # Sparser estimate pays less complexity penalty
  dense <- lapply(X, function(x) solve(cov(x)))
  E <- function(th) sum(abs(th) >= 0.001)
  expect_equal(ggm_bic(X, dense) - dina_vBIC(X, list(theta = dense)), 0)
  expect_equal(ggm_bic(list(X[[1]]), list(diag(4))),
               sum(diag(cov(X[[1]]))) + log(25) * 4 / 25)
  expect_equal(E(diag(4)), 4)                                     # diagonal counted
})

test_that("ggm_bic returns NA for invalid precision matrices (changed from DiNA)", {
  set.seed(1)
  X <- list(matrix(rnorm(40), 10, 4))
  bad <- list(diag(c(1, 1, 1, -1)))
  expect_true(is.na(ggm_bic(X, bad)))
  expect_false(is.na(dina_vBIC(X, list(theta = bad))))           # DiNA ignored the sign
  expect_true(is.na(ggm_bic(X, list(matrix(NA_real_, 4, 4)))))
  Xna <- X; Xna[[1]][1, 1] <- NA
  expect_true(is.na(ggm_bic(Xna, list(diag(4)))))
  expect_error(ggm_bic(X, list(diag(4), diag(4))), "same, non-zero length")
})
