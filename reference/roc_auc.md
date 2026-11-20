# ROC curve and AUC for edge recovery

Computes the ROC curve of a score matrix against a true adjacency
matrix, using the upper triangle. The curve steps only at distinct score
values, so tied scores are handled correctly: a block of tied edges
gives one diagonal segment, i.e. each tied positive/negative pair counts
1/2. The AUC therefore equals the Mann-Whitney statistic and does not
depend on the order of the edges.

## Usage

``` r
roc_auc(score_mat, truth_mat)
```

## Arguments

- score_mat:

  Numeric p x p matrix of edge scores (higher = more likely an edge), or
  a vector of edge scores.

- truth_mat:

  p x p matrix of true edges (non-zero = edge), or a vector matching
  `score_mat`.

## Value

A list with `FPR` and `TPR` (the curve, from (0, 0) to (1, 1)),
`thresholds` (the score at each point; `Inf` for the origin) and `auc`.
If the truth has no edges or no non-edges, the curve is the diagonal and
`auc` is `NA`.

## Details

Ties matter for sparse estimators, which give many edges a score of
exactly 0 (and PIPs that underflow to exactly 0 or 1). An implementation
that accumulates one edge at a time in sorted order breaks ties by
matrix position and can be off by several hundredths of AUC.

Edges with a missing score or truth are dropped.

## Examples

``` r
truth <- matrix(0, 5, 5); truth[1, 2] <- truth[2, 3] <- 1
truth <- truth + t(truth)
score <- matrix(0, 5, 5); score[1, 2] <- 0.9; score[3, 4] <- 0.9
score <- score + t(score)
roc_auc(score, truth)$auc   # ties at 0.9 and at 0 count 1/2
#> [1] 0.6875
```
