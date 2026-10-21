#' Threshold a matrix to an adjacency matrix
#'
#' Declares an edge where \eqn{|m_{ij}| > threshold}, e.g. to turn an
#' estimated precision or partial-correlation matrix into a graph. Applied to
#' a difference \code{theta1 - theta2}, it gives magnitude-based
#' differential edges.
#'
#' @param mat A square numeric matrix.
#' @param threshold Non-negative number. Default 0.01.
#'
#' @return An integer 0/1 matrix with zero diagonal and the dimnames of
#'   \code{mat}. Missing entries stay \code{NA}.
#' @seealso [diff_adjacency()], [confusion_at_threshold()]
#' @export
#'
#' @examples
#' Omega <- matrix(c(2, -0.8, 0.005, -0.8, 2, -0.5, 0.005, -0.5, 1.5), 3, 3)
#' threshold_to_adjacency(Omega)
threshold_to_adjacency <- function(mat, threshold = 0.01) {
  mat <- as.matrix(mat)
  if (nrow(mat) != ncol(mat)) stop("mat must be square")
  adj <- (abs(mat) > threshold) * 1L
  diag(adj) <- 0L
  adj
}

#' Structural differential edges between two graphs
#'
#' An edge is differential when it is present in exactly one of the two
#' graphs (exclusive or). For differential edges defined by the size of the
#' change instead, use \code{threshold_to_adjacency(theta1 - theta2,
#' threshold)}.
#'
#' @param adj1,adj2 Square matrices of the same dimensions; non-zero entries
#'   are edges.
#'
#' @return An integer 0/1 matrix with zero diagonal and the dimnames of
#'   \code{adj1}. Missing entries stay \code{NA}.
#' @seealso [threshold_to_adjacency()]
#' @export
#'
#' @examples
#' a1 <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3)
#' a2 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
#' diff_adjacency(a1, a2)
diff_adjacency <- function(adj1, adj2) {
  adj1 <- as.matrix(adj1); adj2 <- as.matrix(adj2)
  if (nrow(adj1) != ncol(adj1)) stop("adj1 must be square")
  if (!identical(dim(adj1), dim(adj2))) stop("adj1 and adj2 must have the same dimensions")
  out <- xor(adj1 != 0, adj2 != 0) * 1L
  dim(out) <- dim(adj1)
  dimnames(out) <- dimnames(adj1)
  diag(out) <- 0L
  out
}

#' Symmetric matrix from an upper triangle
#'
#' Inverse of \code{m[upper.tri(m)]}: rebuilds a symmetric p x p matrix from
#' its upper triangle in \code{upper.tri()} order (the order of
#' [edge_index()]). Useful for methods that return one value per edge, e.g.
#' a vector of partial correlations.
#'
#' @param x Numeric vector of length \eqn{p(p - 1)/2}.
#' @param varnames Optional character vector of length \code{p} used as
#'   dimnames.
#' @param diag Value for the diagonal. Default 1 (partial correlations).
#'
#' @return A symmetric p x p matrix.
#' @seealso [edge_index()]
#' @export
#'
#' @examples
#' upper_tri_to_matrix(c(0.5, 0.1, -0.3), c("a", "b", "c"))
upper_tri_to_matrix <- function(x, varnames = NULL, diag = 1) {
  p <- (1 + sqrt(1 + 8 * length(x))) / 2
  if (length(x) < 1 || abs(p - round(p)) > 1e-8)
    stop("length(x) must be p * (p - 1) / 2 for an integer p >= 2")
  p <- as.integer(round(p))
  if (!is.null(varnames) && length(varnames) != p) stop("varnames must have length p")
  m <- matrix(0, p, p)
  m[upper.tri(m)] <- x
  m <- m + t(m)
  diag(m) <- diag
  if (!is.null(varnames)) dimnames(m) <- list(varnames, varnames)
  m
}

#' Table of differential edges between two estimates
#'
#' Lists the edges whose estimates differ by more than \code{threshold}
#' between two groups, in the shared edge-table format.
#'
#' @param mat1,mat2 Symmetric p x p matrices of edge estimates for groups 1
#'   and 2 (e.g. partial correlations).
#' @param threshold Non-negative number; edges with
#'   \eqn{|mat1 - mat2| > threshold} are kept. Default 0.05. Use \code{-Inf}
#'   to keep every edge.
#' @param varnames Optional variable names; defaults to the row names of
#'   \code{mat1}, then \code{V1, ..., Vp}.
#'
#' @return A data frame with columns \code{i}, \code{j}, \code{node1},
#'   \code{node2} (as in [edge_index()]), \code{est_g1}, \code{est_g2} and
#'   \code{diff} (\code{est_g1 - est_g2}), in \code{upper.tri()} order. Edges
#'   with a missing difference are dropped. Sort by \code{-abs(diff)} to list
#'   the largest changes first.
#' @seealso [edge_index()], [pooled_perm_pvalues()]
#' @export
#'
#' @examples
#' P1 <- matrix(c(1, 0.4, 0, 0.4, 1, 0.3, 0, 0.3, 1), 3, 3)
#' P2 <- matrix(c(1, 0.1, 0, 0.1, 1, 0.28, 0, 0.28, 1), 3, 3)
#' diff_edge_table(P1, P2)
diff_edge_table <- function(mat1, mat2, threshold = 0.05, varnames = NULL) {
  mat1 <- as.matrix(mat1); mat2 <- as.matrix(mat2)
  if (nrow(mat1) != ncol(mat1)) stop("mat1 must be square")
  if (!identical(dim(mat1), dim(mat2))) stop("mat1 and mat2 must have the same dimensions")
  ut <- upper.tri(mat1)
  out <- edge_index(nrow(mat1), varnames %||% rownames(mat1))
  out$est_g1 <- mat1[ut]
  out$est_g2 <- mat2[ut]
  out$diff <- out$est_g1 - out$est_g2
  keep <- !is.na(out$diff) & abs(out$diff) > threshold
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}
