# ggmtools

Shared tools for R packages that estimate multiple Gaussian graphical
models (GGMs), currently

- [spikeyglass](https://github.com/mljaniczek/spikeyglass): Bayesian
  joint spike-and-slab graphical lasso (Li, McCormick & Clark 2019),
  fitted by EM;
- [multiGGM](https://github.com/mljaniczek/multiGGMr): Bayesian multiple
  GGMs by MCMC: G-Wishart with a Markov random field prior (Peterson,
  Stingo & Vannucci 2015), scalable spike-and-slab (scale-multiGGM;
  Shaddox et al. 2018), multi-platform (Shaddox et al. 2020) and the
  multiple graphical horseshoe (mGHS; Busatto & Stingo 2025).

Both packages import ggmtools and re-export what users call, so the
functions are the same objects whichever package is attached and one
downstream script works on fits from either.

Most functions take plain matrices (precision, partial-correlation,
score or adjacency matrices) or vectors of p-values and probabilities,
not fitted objects. So they work just as well with any other estimator,
such as the joint graphical lasso (JGL), the graphical lasso or a
random-forest importance score: see [Using ggmtools with other
methods](https://mljaniczek.github.io/ggmtools/articles/other-methods.html).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("mljaniczek/ggmtools")
```

## What is in it

| What | Functions |
|----|----|
| Shared edge-table format | [`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md), [`as_pip_table()`](https://mljaniczek.github.io/ggmtools/reference/as_pip_table.md), [`diff_edge_table()`](https://mljaniczek.github.io/ggmtools/reference/diff_edge_table.md) |
| Generics with methods in each package | [`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md) (long table of PIPs), [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md) (one partial-correlation matrix per group; ggmtools also provides the method for JGL fits) |
| Evaluation against a known graph | [`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md), [`pr_auc()`](https://mljaniczek.github.io/ggmtools/reference/pr_auc.md) (tie-correct ROC and PR curves), [`plot_roc()`](https://mljaniczek.github.io/ggmtools/reference/plot_roc.md), [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md) (TPR, FPR, precision, F1, FDR, MCC), [`topk_overlap()`](https://mljaniczek.github.io/ggmtools/reference/topk_overlap.md) |
| Estimation error and tuning | [`relative_frobenius()`](https://mljaniczek.github.io/ggmtools/reference/relative_frobenius.md), [`kl_divergence()`](https://mljaniczek.github.io/ggmtools/reference/kl_divergence.md), [`ggm_bic()`](https://mljaniczek.github.io/ggmtools/reference/ggm_bic.md) |
| Edge selection | [`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md) (Bayesian FDR from posterior probabilities), [`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md) (pooled permutation null) |
| Matrix utilities | [`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md), [`threshold_to_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/threshold_to_adjacency.md), [`diff_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/diff_adjacency.md), [`upper_tri_to_matrix()`](https://mljaniczek.github.io/ggmtools/reference/upper_tri_to_matrix.md) |
| Data | `brca_pi3k`: TCGA BRCA PI3K-pathway expression, Normal vs Tumor (100 x 28 each) |

## Example: differential edges in the TCGA data

`brca_pi3k` holds standardized expression of 28 PI3K-pathway genes in
100 normal and 100 tumour breast samples. With n \> p, the inverse
sample covariance gives a simple (unregularized) estimate of each
group’s partial correlations.

``` r

library(ggmtools)

pcor <- lapply(brca_pi3k, function(x) precision_to_pcor(solve(cov(x))))
d <- diff_edge_table(pcor$Normal, pcor$Tumor, threshold = 0.3)
head(d[order(-abs(d$diff)), ])
#>     i  j  node1   node2     est_g1      est_g2       diff
#> 37  2 21 AKT1S1   NPRL3  0.4344316 -0.16654190  0.6009735
#> 33  8 18   MTOR  DEPDC5 -0.1053410  0.47261097 -0.5779520
#> 23  3 15   AKT2    RPS6 -0.5468767  0.02156526 -0.5684420
#> 51 13 24   RHEB PPP2R1A  0.3695557 -0.18810932  0.5576650
#> 47  4 24   AKT3 PPP2R1A  0.3664210 -0.16630078  0.5327218
#> 5   3  6   AKT2 MAPKAP1  0.3845418 -0.14658324  0.5311250
```

Which of these differences are larger than chance? Permuting the group
labels gives a null distribution of \|difference\|; pooling it over all
edges gives fine-grained p-values.

``` r

set.seed(1)
X <- rbind(brca_pi3k$Normal, brca_pi3k$Tumor)
diff_score <- function(rows) {
  p1 <- precision_to_pcor(solve(cov(X[rows, ])))
  p2 <- precision_to_pcor(solve(cov(X[-rows, ])))
  abs(p1 - p2)
}
observed <- diff_score(1:100)
null <- replicate(100, diff_score(sample(200, 100)), simplify = FALSE)

pp <- pooled_perm_pvalues(observed, null)
edges <- cbind(diff_edge_table(pcor$Normal, pcor$Tumor, threshold = -Inf),
               pvalue = pp$pvalue, p_adj = pp$p_adj)
head(edges[order(edges$pvalue), c("node1", "node2", "diff", "pvalue", "p_adj")])
#>      node1   node2       diff      pvalue    p_adj
#> 192 AKT1S1   NPRL3  0.6009735 0.004550144 0.604984
#> 144   MTOR  DEPDC5 -0.5779520 0.006507764 0.604984
#> 94    AKT2    RPS6 -0.5684420 0.007380757 0.604984
#> 266   RHEB PPP2R1A  0.5576650 0.008491839 0.604984
#> 257   AKT3 PPP2R1A  0.5327218 0.011851538 0.604984
#> 13    AKT2 MAPKAP1  0.5311250 0.012010264 0.604984
sum(edges$p_adj < 0.05)
#> [1] 0
```

With 100 samples per group and no regularization, no difference survives
the BH correction over all 378 edges. The joint estimators in
spikeyglass and multiGGM borrow strength across groups and are the tools
to use for a real analysis; the same two calls then test their
estimates.

A Bayesian method gives posterior probabilities instead, e.g. P(edge
differs).
[`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md)
selects the largest set of edges whose posterior expected FDR is below a
target (tied probabilities are kept together):

``` r

prob <- c(0.999, 0.99, 0.97, 0.97, 0.80, 0.40, 0.05)
bayes_fdr(prob, alpha = 0.05)
#> [1]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE
#> attr(,"expected_fdr")
#> [1] 0.01775
#> attr(,"n_selected")
#> [1] 4
#> attr(,"cutoff")
#> [1] 0.97
#> attr(,"alpha")
#> [1] 0.05
```

## Example: evaluating an estimate against a known graph

In a simulation the true precision matrix is known. Here, a chain graph
on 20 variables, estimated from only 40 samples:

``` r

set.seed(3)
p <- 20
Omega <- diag(p)
Omega[cbind(1:(p - 1), 2:p)] <- Omega[cbind(2:p, 1:(p - 1))] <- 0.3
X <- matrix(rnorm(40 * p), 40) %*% chol(solve(Omega))

Omega_hat <- solve(cov(X))
score <- abs(precision_to_pcor(Omega_hat))   # edge score: |partial correlation|
truth <- threshold_to_adjacency(Omega, threshold = 0)

roc_auc(score, truth)$auc
#> [1] 0.8024007
pr_auc(score, truth)$auc
#> [1] 0.3226991
cm <- confusion_at_threshold(score, truth, threshold = 0.3)
unlist(cm[c("TP", "FP", "TPR", "FDR", "MCC")])
#>         TP         FP        TPR        FDR        MCC 
#> 11.0000000 23.0000000  0.5789474  0.6764706  0.3478489
kl_divergence(Omega_hat, solve(Omega))
#> [1] 6.767717
relative_frobenius(Omega_hat, Omega)
#> [1] 2.173401
```

ROC and PR curves step only at distinct score values, so sparse
estimators, which tie many edges at exactly 0, get an AUC that does not
depend on the order of the variables. For a sparse estimate, also
compare `pr_auc(..., method = "average_precision")`: the default
trapezoid draws a straight line across a large tied block, which is
optimistic.

``` r

sparse <- score * (score > 0.3)      # e.g. a thresholded or lasso-type estimate
mean(sparse[upper.tri(sparse)] == 0) # share of edges tied at 0
#> [1] 0.8210526
c(trapezoid = pr_auc(sparse, truth)$auc,
  average_precision = pr_auc(sparse, truth, method = "average_precision")$auc)
#>         trapezoid average_precision 
#>         0.3188026         0.2854237
```

## Edge tables

Every edge table has the columns `i`, `j`, `node1`, `node2` (one row per
pair i \< j, in [`upper.tri()`](https://rdrr.io/r/base/lower.tri.html)
order), followed by method-specific columns.
[`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md)
adds `group` and `pip`. Differential-edge tables
([`diff_edge_table()`](https://mljaniczek.github.io/ggmtools/reference/diff_edge_table.md),
`spikeyglass::diff_edges()`) add `est_g1`, `est_g2`, `diff` and
inference columns, using `lower`/`upper` for intervals.

## References

- Busatto, C. & Stingo, F.C. (2025). Inference of multiple
  high-dimensional networks with the Graphical Horseshoe prior. *Journal
  of Computational and Graphical Statistics*.
  <doi:10.1080/10618600.2024.2435480>
- Danaher, P., Wang, P. & Witten, D.M. (2014). The joint graphical lasso
  for inverse covariance estimation across multiple classes. *Journal of
  the Royal Statistical Society: Series B*, 76(2), 373-397.
- Davis, J. & Goadrich, M. (2006). The relationship between
  precision-recall and ROC curves. *Proceedings of the 23rd
  International Conference on Machine Learning*, 233-240.
- Li, Z.R., McCormick, T.H. & Clark, S.J. (2019). Bayesian joint
  spike-and-slab graphical lasso. *Proceedings of Machine Learning
  Research*, 97, 3877-3885.
- Newton, M.A., Noueiry, A., Sarkar, D. & Ahlquist, P. (2004). Detecting
  differential gene expression with a semiparametric hierarchical
  mixture method. *Biostatistics*, 5(2), 155-176.
- Peterson, C.B., Stingo, F.C. & Vannucci, M. (2015). Bayesian inference
  of multiple Gaussian graphical models. *Journal of the American
  Statistical Association*, 110(509), 159-174.
- Shaddox, E., Stingo, F.C., Peterson, C.B., et al. (2018). A Bayesian
  approach for learning gene networks underlying disease severity in
  COPD. *Statistics in Biosciences*, 10(1), 59-85.
- Shaddox, E., Peterson, C.B., Stingo, F.C., et al. (2020). Bayesian
  inference of networks across multiple sample groups and data types.
  *Biostatistics*, 21(3), 561-576.
