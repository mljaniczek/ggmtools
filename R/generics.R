#' Posterior inclusion probabilities as a long edge table
#'
#' Generic for extracting edge-level posterior inclusion probabilities (PIPs)
#' from a fitted multiple-GGM model, in a format shared by all packages that
#' implement it, so that one downstream script works on any of their fits.
#'
#' @param fit A fitted model object.
#' @param ... Passed to methods (e.g. which ladder step or chain to use).
#'
#' @return A data frame with one row per edge \eqn{i < j} and group, ordered by
#'   group and then in \code{upper.tri()} order, with columns
#'   \describe{
#'     \item{i, j}{Variable indices, \eqn{i < j}.}
#'     \item{node1, node2}{Variable names.}
#'     \item{group}{Group name (character).}
#'     \item{pip}{Posterior probability that the edge is in that group's
#'       graph.}
#'   }
#'   Methods document how their PIPs are defined. For example, in
#'   \pkg{spikeyglass} the edge indicator is shared by all groups, so
#'   \code{pip} is identical across groups; in \pkg{multiGGMr} each group has
#'   its own graph.
#'
#' @seealso [pcor_estimate()], [edge_index()]
#' @export
pip_table <- function(fit, ...) UseMethod("pip_table")

#' Point estimates of partial correlations
#'
#' Generic for extracting one partial-correlation matrix per group from a
#' fitted multiple-GGM model.
#'
#' @param fit A fitted model object.
#' @param ... Passed to methods.
#'
#' @return A named list of K symmetric p x p matrices with unit diagonal and
#'   variable names as dimnames. The names are the group names, or
#'   \code{"Group_1"}, \code{"Group_2"}, ... if the groups are unnamed. Methods
#'   document which estimate they return. For example, \pkg{spikeyglass}
#'   returns the MAP estimate and \pkg{multiGGMr} the posterior mean.
#'
#' @seealso [pip_table()], [precision_to_pcor()]
#' @export
pcor_estimate <- function(fit, ...) UseMethod("pcor_estimate")

#' Edge index in the shared edge-table format
#'
#' The first four columns of every edge table in the shared format: one row
#' per variable pair \eqn{i < j}, in \code{upper.tri()} order (column-major,
#' i.e. (1, 2), (1, 3), (2, 3), (1, 4), ...).
#'
#' @param p Number of variables.
#' @param varnames Optional character vector of length \code{p}. Defaults to
#'   \code{V1, ..., Vp}.
#'
#' @return A data frame with columns \code{i}, \code{j}, \code{node1},
#'   \code{node2}.
#' @export
#'
#' @examples
#' edge_index(4, c("a", "b", "c", "d"))
edge_index <- function(p, varnames = NULL) {
  p <- as.integer(p)
  if (length(p) != 1 || is.na(p) || p < 2) stop("p must be an integer >= 2")
  if (is.null(varnames)) varnames <- paste0("V", seq_len(p))
  if (length(varnames) != p) stop("varnames must have length p")
  ij <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  data.frame(i = ij[, 1], j = ij[, 2],
             node1 = as.character(varnames[ij[, 1]]),
             node2 = as.character(varnames[ij[, 2]]),
             stringsAsFactors = FALSE)
}

#' Build a long PIP table from per-group PIP matrices
#'
#' Helper for implementing [pip_table()] methods.
#'
#' @param pip_list A list of K symmetric p x p PIP matrices (one per group;
#'   repeat the same matrix if the PIP is shared across groups).
#' @param varnames Optional variable names; defaults to the dimnames of the
#'   first matrix, then \code{V1, ..., Vp}.
#' @param groupnames Optional group names; defaults to \code{names(pip_list)},
#'   then \code{Group_1, ..., Group_K}.
#'
#' @return A data frame in the format described in [pip_table()].
#' @export
#'
#' @examples
#' m <- matrix(c(0, .9, .1, .9, 0, .4, .1, .4, 0), 3, 3)
#' as_pip_table(list(A = m, B = m))
as_pip_table <- function(pip_list, varnames = NULL, groupnames = NULL) {
  if (!is.list(pip_list) || length(pip_list) < 1) stop("pip_list must be a non-empty list")
  p <- nrow(pip_list[[1]])
  if (is.null(varnames)) varnames <- rownames(pip_list[[1]])
  groupnames <- default_group_names(groupnames %||% names(pip_list), length(pip_list))
  idx <- edge_index(p, varnames)
  ut <- upper.tri(pip_list[[1]])
  out <- do.call(rbind, lapply(seq_along(pip_list), function(k) {
    data.frame(idx, group = groupnames[k], pip = unname(pip_list[[k]][ut]),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

default_group_names <- function(groupnames, K) {
  if (is.null(groupnames) || length(groupnames) != K || any(groupnames == ""))
    groupnames <- paste0("Group_", seq_len(K))
  as.character(groupnames)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
