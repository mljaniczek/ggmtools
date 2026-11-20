# Empirical p-values against a pooled permutation null

Computes one-sided empirical p-values for observed edge scores (larger =
more extreme, e.g. \\\|\hat\theta^{(1)}\_{ij} -
\hat\theta^{(2)}\_{ij}\|\\) against a null distribution pooled over all
edges and permutations, and Benjamini-Hochberg adjusts them.

## Usage

``` r
pooled_perm_pvalues(observed, null, tol = sqrt(.Machine$double.eps))
```

## Arguments

- observed:

  Numeric vector of observed scores, or a square matrix (the upper
  triangle is used, in
  [`upper.tri()`](https://rdrr.io/r/base/lower.tri.html) order as in
  [`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md)).

- null:

  Null scores: a numeric vector or matrix, or a list of them. All finite
  values are pooled. If `observed` is a matrix, list elements with the
  same dimensions contribute their upper triangle.

- tol:

  Relative tolerance: a null value counts as at least as large as \\x\\
  if it is \\\> x - tol \cdot \max(1, \|x\|)\\, so that ties up to
  rounding error count as exceedances.

## Value

A list with `pvalue`, \\(1 + \\\\null \ge obs\\) / (1 + n\_{null})\\,
and `p_adj` (BH-adjusted), one per observed edge. Missing observed
scores give `NA`. If the pool has no finite values, all p-values are
`NA` with a warning.

## Details

Pooling gives a p-value resolution of about \\1 / (1 + B m)\\ for \\B\\
permutations and \\m\\ edges, instead of \\1 / (1 + B)\\ for per-edge
nulls, which is fine enough to survive a BH correction over many edges.
It assumes the per-edge null distributions are approximately
exchangeable, which is reasonable for label-permutation nulls of
comparable (e.g. standardized) scores.

## See also

[`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md)

## Examples

``` r
set.seed(1)
p <- 6
null <- replicate(50, abs(matrix(rnorm(p * p), p)), simplify = FALSE)
obs <- abs(matrix(rnorm(p * p), p)); obs[1, 2] <- 4
pp <- pooled_perm_pvalues(obs, null)
cbind(edge_index(p), pvalue = pp$pvalue, p_adj = pp$p_adj)[1:3, ]
#>   i j node1 node2      pvalue      p_adj
#> 1 1 2    V1    V2 0.001331558 0.01997337
#> 2 1 3    V1    V3 0.924101198 0.96404794
#> 3 2 3    V2    V3 0.382157124 0.84458817
```
