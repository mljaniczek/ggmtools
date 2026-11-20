# Confusion metrics for edge selection at a threshold

Compares a score matrix (e.g. posterior inclusion probabilities) with a
true adjacency matrix. An edge is predicted when its score is
`>= threshold`. A 0/1 (or logical) estimated adjacency matrix can be
passed as `score_mat` with the default threshold.

## Usage

``` r
confusion_at_threshold(
  score_mat,
  truth_mat,
  threshold = 0.5,
  upper_only = TRUE
)
```

## Arguments

- score_mat:

  Numeric p x p matrix of edge scores, or a vector of edge scores (e.g.
  an upper triangle) of the same length as `truth_mat`.

- truth_mat:

  p x p matrix of true edges (non-zero = edge), or a vector matching
  `score_mat`.

- threshold:

  Numeric threshold. Default 0.5.

- upper_only:

  Logical. If `TRUE` (default), only the upper triangle is used (each
  undirected edge once); otherwise all off-diagonal entries. Ignored for
  vectors.

## Value

A list with `TP`, `FP`, `TN`, `FN`, `TPR` (sensitivity, power), `FPR`,
`precision`, `F1`, `FDR` (`1 - precision`) and `MCC` (Matthews
correlation coefficient). Rates whose denominator is zero are `NA`; in
particular `FDR` is `NA` when no edge is predicted and `MCC` is `NA`
when any margin of the confusion table is 0.

## Details

Edges with a missing score or truth are dropped.

## Examples

``` r
truth <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
score <- matrix(c(0, .9, .6, .9, 0, .2, .6, .2, 0), 3, 3)
confusion_at_threshold(score, truth, threshold = 0.5)
#> $TP
#> [1] 1
#> 
#> $FP
#> [1] 1
#> 
#> $TN
#> [1] 0
#> 
#> $FN
#> [1] 1
#> 
#> $TPR
#> [1] 0.5
#> 
#> $FPR
#> [1] 1
#> 
#> $precision
#> [1] 0.5
#> 
#> $F1
#> [1] 0.5
#> 
#> $FDR
#> [1] 0.5
#> 
#> $MCC
#> [1] -0.5
#> 

# A 0/1 estimated adjacency matrix
confusion_at_threshold(threshold_to_adjacency(score, 0.5), truth)$MCC
#> [1] -0.5
```
