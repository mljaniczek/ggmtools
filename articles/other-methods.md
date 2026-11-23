# Using ggmtools with other methods

ggmtools has two layers:

- **Generics**
  ([`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md),
  [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md))
  that a package implements for its own fitted objects. spikeyglass and
  multiGGM provide methods, and ggmtools itself provides
  [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md)
  for JGL fits, so the same downstream code runs on fits from any of
  them.
- **Matrix-level functions** for everything else: evaluation, losses,
  tuning, edge selection and edge tables. These take precision,
  partial-correlation, score or adjacency matrices, or vectors of
  probabilities and p-values, never a fitted object. They work with any
  estimator once its output is pulled out as a matrix.

This article shows the second layer with the joint graphical lasso (JGL;
Danaher, Wang & Witten 2014), then how to give a new method the first
layer.

## Which functions apply to which output

| The method returns | Examples | ggmtools functions |
|----|----|----|
| Precision matrices, one per group | JGL (fused or group; [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md) works on its fits), graphical lasso per group, spikeyglass, inverse covariance | [`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md), [`threshold_to_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/threshold_to_adjacency.md), [`kl_divergence()`](https://mljaniczek.github.io/ggmtools/reference/kl_divergence.md), [`relative_frobenius()`](https://mljaniczek.github.io/ggmtools/reference/relative_frobenius.md), [`ggm_bic()`](https://mljaniczek.github.io/ggmtools/reference/ggm_bic.md) (tuning), then everything below |
| Partial correlations, as matrices or upper-triangle vectors | multiGGM [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md), DiNGO | [`upper_tri_to_matrix()`](https://mljaniczek.github.io/ggmtools/reference/upper_tri_to_matrix.md), [`diff_edge_table()`](https://mljaniczek.github.io/ggmtools/reference/diff_edge_table.md), [`threshold_to_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/threshold_to_adjacency.md) |
| Edge scores (higher = more likely an edge or a difference) | \|partial correlation\|, random-forest importance (JRF), D-trace \|Δ\|, test statistics | [`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md), [`pr_auc()`](https://mljaniczek.github.io/ggmtools/reference/pr_auc.md), [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md), [`topk_overlap()`](https://mljaniczek.github.io/ggmtools/reference/topk_overlap.md) |
| Posterior edge probabilities | multiGGM (MRF priors and mGHS), spikeyglass [`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md) | [`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md), and the score functions above |
| Permutation or bootstrap null scores | any method refitted on permuted labels | [`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md) |
| Estimated graphs (0/1) | any thresholded or selected graph | [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md), [`diff_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/diff_adjacency.md) |

## Example: the joint graphical lasso

The JGL package was archived on CRAN in July 2025; install it from the
CRAN archive or with `remotes::install_github("cran/JGL")`. It returns a
list whose `theta` element holds one precision matrix per group. Tune
the two penalties with
[`ggm_bic()`](https://mljaniczek.github.io/ggmtools/reference/ggm_bic.md),
the BIC used with JGL in the Koyejo lab code (smaller is better):

``` r

library(ggmtools)
library(JGL)

grid <- expand.grid(lambda1 = c(0.05, 0.1, 0.2, 0.3),
                    lambda2 = c(0.01, 0.05, 0.1))
grid$bic <- apply(grid, 1, function(g) {
  fit <- JGL(brca_pi3k, penalty = "fused", lambda1 = g[["lambda1"]],
             lambda2 = g[["lambda2"]], return.whole.theta = TRUE)
  ggm_bic(brca_pi3k, fit$theta)
})
best <- grid[which.min(grid$bic), ]
best
#>   lambda1 lambda2      bic
#> 2     0.1    0.01 29.43266
```

Refit at the chosen penalties. ggmtools has a
[`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md)
method for JGL fits, which converts each precision matrix with
[`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md).
JGL does not store group names, so pass them, and fit with
`return.whole.theta = TRUE` so that every variable is kept:

``` r

fit <- JGL(brca_pi3k, penalty = "fused", lambda1 = best$lambda1,
           lambda2 = best$lambda2, return.whole.theta = TRUE)
genes <- colnames(brca_pi3k$Normal)
pcor <- pcor_estimate(fit, groupnames = names(brca_pi3k))

# Number of edges in each group's graph
sapply(fit$theta, function(theta) sum(threshold_to_adjacency(theta, 1e-8)) / 2)
#> [1] 123 149
```

The fused penalty sets many differences to exactly 0.
[`diff_edge_table()`](https://mljaniczek.github.io/ggmtools/reference/diff_edge_table.md)
lists the rest in the shared edge-table format:

``` r

d <- diff_edge_table(pcor$Normal, pcor$Tumor, threshold = 1e-8, varnames = genes)
nrow(d)
#> [1] 201
head(d[order(-abs(d$diff)), ])
#>      i  j  node1   node2      est_g1     est_g2       diff
#> 7    3  6   AKT2 MAPKAP1  0.38552727  0.0000000  0.3855273
#> 90   9 19   PDK1  INPP4B -0.04251113 -0.3683474  0.3258362
#> 80   8 18   MTOR  DEPDC5  0.02470883  0.3490456 -0.3243367
#> 107  2 21 AKT1S1   NPRL3  0.29674004  0.0000000  0.2967400
#> 169 17 26  RPTOR   STK11  0.29551205  0.0000000  0.2955121
#> 180 10 27 PIK3CA    TSC1  0.27580379  0.0000000  0.2758038
```

Do JGL and the unregularized estimate agree on which edges change most?
JGL ties every edge it fuses at a difference of exactly 0, so
[`topk_overlap()`](https://mljaniczek.github.io/ggmtools/reference/topk_overlap.md)
counts those ties fairly instead of by edge position:

``` r

pcor_inv <- lapply(brca_pi3k, function(x) precision_to_pcor(solve(cov(x))))
k <- c(10, 25, 50, 100, 200)
data.frame(k = k,
           overlap = topk_overlap(abs(pcor$Normal - pcor$Tumor),
                                  abs(pcor_inv$Normal - pcor_inv$Tumor), k))
#>     k overlap
#> 1  10       3
#> 2  25       4
#> 3  50      15
#> 4 100      35
#> 5 200     110
```

In a simulation, where the true graphs are known, the same estimates go
straight into
[`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md),
[`pr_auc()`](https://mljaniczek.github.io/ggmtools/reference/pr_auc.md)
and
[`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md),
with `abs(pcor$Normal)` (or
`threshold_to_adjacency(fit$theta[[1]], 1e-8)`) as the score and the
true adjacency matrix as the truth.

## Example: the multiple graphical horseshoe (mGHS)

mGHS (Busatto & Stingo 2025) is fitted in multiGGM with
`multiggm_mcmc(method = "mghs")`. The fit is a multiGGM fit, so the
ggmtools generics work directly. The horseshoe prior has no point mass
at zero, so multiGGM defines the edge probability from the shrinkage
weights, $`P(w_{k,ij} > 0.5)`$ with $`w = \lambda^2 / (1 + \lambda^2)`$;
see the multiGGM documentation.

``` r

library(multiGGM)
fit <- multiggm_mcmc(data_list = brca_pi3k, method = "mghs")

pips <- pip_table(fit)       # long table: i, j, node1, node2, group, pip
pcor <- pcor_estimate(fit)   # list of partial-correlation matrices

normal <- pips[pips$group == "Normal", ]
normal$selected <- bayes_fdr(normal$pip, alpha = 0.05)
```

## Giving a new method the generics

[`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md)
and
[`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md)
form a common interface. Supporting a new method means writing an
*adapter*: S3 methods that convert its fitted object into the shared
output. ggmtools ships the adapter for JGL fits (class `"jgl"`), and it
is a template for others. Here it is for a hypothetical method whose
fits have class `"mymethod"` and store a list of precision matrices in
`fit$precision`:

``` r

pcor_estimate.mymethod <- function(fit, varnames = NULL, groupnames = NULL, ...) {
  out <- lapply(fit$precision, precision_to_pcor)
  p <- nrow(out[[1]])
  varnames <- if (!is.null(varnames)) varnames else
    if (!is.null(rownames(out[[1]]))) rownames(out[[1]]) else paste0("V", seq_len(p))
  out <- lapply(out, function(m) { dimnames(m) <- list(varnames, varnames); m })
  names(out) <- if (!is.null(groupnames)) groupnames else
    paste0("Group_", seq_along(out))
  out
}
```

In a script, defining the function is enough for dispatch. In a package,
register it with `S3method(ggmtools::pcor_estimate, mymethod)` in the
NAMESPACE (with roxygen: `@exportS3Method ggmtools::pcor_estimate`). A
method that returns posterior probabilities can build its
[`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md)
output with
[`as_pip_table()`](https://mljaniczek.github.io/ggmtools/reference/as_pip_table.md).

## References

- Busatto, C. & Stingo, F.C. (2025). Inference of multiple
  high-dimensional networks with the Graphical Horseshoe prior. *Journal
  of Computational and Graphical Statistics*.
  [doi:10.1080/10618600.2024.2435480](https://doi.org/10.1080/10618600.2024.2435480)
- Danaher, P., Wang, P. & Witten, D.M. (2014). The joint graphical lasso
  for inverse covariance estimation across multiple classes. *Journal of
  the Royal Statistical Society: Series B*, 76(2), 373-397.
