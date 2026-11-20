# Precision-recall curve and AUPR for edge recovery

Computes the precision-recall (PR) curve of a score matrix against a
true adjacency matrix, using the upper triangle. As in
[`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md),
the curve steps only at distinct score values (one point per block of
tied scores), so the result does not depend on the order of the edges.

## Usage

``` r
pr_auc(
  score_mat,
  truth_mat,
  method = c("trapezoid", "average_precision", "interpolated")
)
```

## Arguments

- score_mat:

  Numeric p x p matrix of edge scores (higher = more likely an edge), or
  a vector of edge scores.

- truth_mat:

  p x p matrix of true edges (non-zero = edge), or a vector matching
  `score_mat`.

- method:

  Summary of the curve: `"trapezoid"` (default), `"average_precision"`
  or `"interpolated"`. See Details.

## Value

A list with `recall` and `precision` (the curve, one point per distinct
score, preceded by the point (0, 1)), `thresholds` (the score at each
point; `Inf` for the first) and `auc`. If the truth has no edges or no
non-edges, `auc` is `NA`.

## Details

Three summaries of the curve are available:

- `"trapezoid"`:

  Trapezoid rule over the curve points, starting from precision 1 at
  recall 0. This matches the tie-corrected `compute_aupr()` of the DiNA
  simulation. A straight line across a large tied block is not
  achievable in PR space, so this is optimistic for sparse estimators
  that tie many edges at 0.

- `"average_precision"`:

  Step-wise average precision, \\\sum_b (R_b - R\_{b-1}) P_b\\ over tie
  blocks \\b\\, as in scikit-learn. A tied block contributes at the
  precision of the whole block.

- `"interpolated"`:

  Exact area under the curve obtained by interpolating true and false
  positives linearly within each tied block (Davis & Goadrich 2006; the
  continuous version, as in the PRROC package). This is the expected PR
  curve under random tie-breaking.

Without ties, `"average_precision"` and `"interpolated"` are close; with
large tied blocks `"trapezoid"` can be much larger than both.

Edges with a missing score or truth are dropped.

## References

Davis, J. & Goadrich, M. (2006). The relationship between
precision-recall and ROC curves. *Proceedings of the 23rd International
Conference on Machine Learning*, 233-240.

## See also

[`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md)

## Examples

``` r
truth <- matrix(0, 5, 5); truth[1, 2] <- truth[2, 3] <- 1
truth <- truth + t(truth)
score <- matrix(0, 5, 5); score[1, 2] <- 0.9; score[3, 4] <- 0.9
score <- score + t(score)
pr_auc(score, truth)$auc
#> [1] 0.55
pr_auc(score, truth, method = "average_precision")$auc
#> [1] 0.35
pr_auc(score, truth, method = "interpolated")$auc
#> [1] 0.3879424
```
