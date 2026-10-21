test_that("threshold_to_adjacency", {
  m <- matrix(c(2, -0.8, 0.005, -0.8, 2, NA, 0.005, NA, 1.5), 3, 3,
              dimnames = list(letters[1:3], letters[1:3]))
  a <- threshold_to_adjacency(m)
  expect_equal(unname(a), matrix(c(0L, 1L, 0L, 1L, 0L, NA, 0L, NA, 0L), 3, 3))
  expect_identical(dimnames(a), dimnames(m))
  expect_type(a, "integer")
  expect_equal(sum(threshold_to_adjacency(m, threshold = 0.8), na.rm = TRUE), 0)
  expect_equal(sum(threshold_to_adjacency(matrix(0, 3, 3))), 0)   # no edges
  expect_error(threshold_to_adjacency(matrix(1, 2, 3)), "square")
})

test_that("diff_adjacency is the XOR of two graphs", {
  a1 <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3)
  a2 <- matrix(c(0, 1, 0, 1, 0, -0.3, 0, -0.3, 0), 3, 3)
  d <- diff_adjacency(a1, a2)
  expect_equal(d, matrix(c(0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L), 3, 3))
  expect_equal(sum(diff_adjacency(a1, a1)), 0)                    # identical graphs
  a3 <- a1; a3[1, 3] <- a3[3, 1] <- NA
  expect_true(is.na(diff_adjacency(a3, a2)[1, 3]))
  expect_error(diff_adjacency(a1, diag(2)), "same dimensions")
})

test_that("upper_tri_to_matrix inverts upper.tri()", {
  m <- upper_tri_to_matrix(c(0.5, 0.1, -0.3), c("a", "b", "c"))
  expect_equal(m["a", "b"], 0.5)
  expect_equal(m["c", "a"], 0.1)
  expect_equal(m["b", "c"], -0.3)
  expect_equal(unname(diag(m)), c(1, 1, 1))
  expect_true(isSymmetric(m))
  x <- runif(10)
  expect_equal(upper_tri_to_matrix(x, diag = 0)[upper.tri(diag(5))], x)
  expect_equal(diag(upper_tri_to_matrix(x, diag = 0)), rep(0, 5))
  expect_true(is.na(upper_tri_to_matrix(c(NA, 1, 2))[2, 1]))
  expect_error(upper_tri_to_matrix(1:4), "p \\* \\(p - 1\\)")
  expect_error(upper_tri_to_matrix(numeric(0)))
  expect_error(upper_tri_to_matrix(1:3, c("a", "b")), "varnames")
})

test_that("diff_edge_table follows the edge-table format", {
  P1 <- matrix(c(1, 0.4, 0, 0.4, 1, 0.3, 0, 0.3, 1), 3, 3,
               dimnames = list(c("x", "y", "z"), c("x", "y", "z")))
  P2 <- matrix(c(1, 0.1, -0.2, 0.1, 1, 0.28, -0.2, 0.28, 1), 3, 3)
  d <- diff_edge_table(P1, P2)
  expect_named(d, c("i", "j", "node1", "node2", "est_g1", "est_g2", "diff"))
  expect_equal(d$i, c(1, 1))
  expect_equal(d$j, c(2, 3))
  expect_equal(d$node1, c("x", "x"))
  expect_equal(d$node2, c("y", "z"))
  expect_equal(d$diff, c(0.3, 0.2))
  expect_equal(nrow(diff_edge_table(P1, P2, threshold = -Inf)), 3)
  expect_equal(nrow(diff_edge_table(P1, P1)), 0)                    # no differences
  expect_named(diff_edge_table(P1, P1), names(d))
  expect_equal(diff_edge_table(unname(P1), P2)$node1, c("V1", "V1"))
  P3 <- P2; P3[1, 2] <- P3[2, 1] <- NA
  expect_equal(diff_edge_table(P1, P3)$j, 3)                         # NA dropped
  expect_error(diff_edge_table(P1, diag(2)), "same dimensions")
})

test_that("precision_to_pcor replaces DiNA's pcor_from_precision", {
  dina_pcor <- function(Theta) {
    d <- sqrt(diag(Theta)); P <- -Theta / outer(d, d); diag(P) <- 1; P
  }
  Om <- matrix(c(2, -0.8, 0, -0.8, 2, -0.5, 0, -0.5, 1.5), 3, 3)
  expect_equal(precision_to_pcor(Om), dina_pcor(Om))
  # An asymmetric estimate (e.g. from neighbourhood selection) is averaged
  Om2 <- Om; Om2[1, 2] <- -0.6
  expect_false(isSymmetric(dina_pcor(Om2)))
  expect_true(isSymmetric(precision_to_pcor(Om2)))
  expect_equal(precision_to_pcor(Om2)[1, 2], (0.6 + 0.8) / 2 / 2)
  # A non-positive diagonal gives NA with a warning instead of NaN
  expect_true(any(is.nan(suppressWarnings(dina_pcor(diag(c(1, -1)))))))
  expect_warning(bad <- precision_to_pcor(diag(c(1, -1))), "non-positive")
  expect_true(all(is.na(bad)) && !any(is.nan(bad)))
})
