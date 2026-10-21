test_that("edge_index follows upper.tri order", {
  e <- edge_index(4, letters[1:4])
  expect_equal(e$i, c(1, 1, 2, 1, 2, 3))
  expect_equal(e$j, c(2, 3, 3, 4, 4, 4))
  expect_equal(e$node1, c("a", "a", "b", "a", "b", "c"))
  expect_equal(edge_index(3)$node2, c("V2", "V3", "V3"))
  expect_error(edge_index(1))
  expect_error(edge_index(3, c("a", "b")))
})

test_that("as_pip_table builds the long format", {
  m1 <- matrix(c(0, .9, .1, .9, 0, .4, .1, .4, 0), 3, 3,
               dimnames = list(c("x", "y", "z"), c("x", "y", "z")))
  m2 <- m1 / 2
  t1 <- as_pip_table(list(m1, m2))
  expect_named(t1, c("i", "j", "node1", "node2", "group", "pip"))
  expect_equal(nrow(t1), 6)
  expect_equal(unique(t1$group), c("Group_1", "Group_2"))
  expect_equal(t1$pip, c(.9, .1, .4, .45, .05, .2))
  expect_equal(t1$node1[1:3], c("x", "x", "y"))
  t2 <- as_pip_table(list(A = m1, B = m2))
  expect_equal(unique(t2$group), c("A", "B"))
})

test_that("generics dispatch to methods", {
  fit <- structure(list(), class = "toyfit")
  pip_table.toyfit <- function(fit, ...) "pip method"
  pcor_estimate.toyfit <- function(fit, ...) "pcor method"
  registerS3method("pip_table", "toyfit", pip_table.toyfit)
  registerS3method("pcor_estimate", "toyfit", pcor_estimate.toyfit)
  expect_equal(pip_table(fit), "pip method")
  expect_equal(pcor_estimate(fit), "pcor method")
})

test_that("precision_to_pcor", {
  Om <- matrix(c(2, -0.8, -0.8, 2), 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  P <- precision_to_pcor(Om)
  expect_equal(P["a", "b"], 0.4)
  expect_equal(unname(diag(P)), c(1, 1))
  expect_warning(bad <- precision_to_pcor(diag(c(1, -1))), "non-positive")
  expect_true(all(is.na(bad)))
})

test_that("brca_pi3k has the documented structure", {
  expect_named(brca_pi3k, c("Normal", "Tumor"))
  expect_equal(dim(brca_pi3k$Normal), c(100, 28))
  expect_equal(dim(brca_pi3k$Tumor), c(100, 28))
  expect_identical(colnames(brca_pi3k$Normal), colnames(brca_pi3k$Tumor))
  expect_false(anyNA(unlist(brca_pi3k)))
  g <- attr(brca_pi3k, "genes")
  expect_equal(g$Gene, colnames(brca_pi3k$Normal))
  expect_true(all(g$pathway == "PI3K"))
})
