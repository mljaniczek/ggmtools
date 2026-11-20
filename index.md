# ggmtools

Shared tools for R packages that estimate multiple Gaussian graphical
models, currently
[spikeyglass](https://github.com/mljaniczek/spikeyglass) (Bayesian
spike-and-slab joint graphical lasso, EM) and
[multiGGMr](https://github.com/mljaniczek/multiGGMr) (Bayesian multiple
GGMs with MRF priors, MCMC).

Both packages import ggmtools and re-export what users call. So the
functions below are the same objects whichever package is attached, and
one downstream script works on fits from either.

| What | Functions |
|----|----|
| Shared edge-table format | [`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md), [`as_pip_table()`](https://mljaniczek.github.io/ggmtools/reference/as_pip_table.md) |
| Generics with methods in each package | [`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md) (long table of PIPs), [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md) (one partial-correlation matrix per group) |
| Evaluation | [`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md) (tie-correct ROC/AUC), [`pr_auc()`](https://mljaniczek.github.io/ggmtools/reference/pr_auc.md) (tie-correct PR/AUPR), [`plot_roc()`](https://mljaniczek.github.io/ggmtools/reference/plot_roc.md), [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md) (incl. FDR, MCC), [`topk_overlap()`](https://mljaniczek.github.io/ggmtools/reference/topk_overlap.md) |
| Estimation error and tuning | [`relative_frobenius()`](https://mljaniczek.github.io/ggmtools/reference/relative_frobenius.md), [`kl_divergence()`](https://mljaniczek.github.io/ggmtools/reference/kl_divergence.md), [`ggm_bic()`](https://mljaniczek.github.io/ggmtools/reference/ggm_bic.md) |
| Inference | [`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md) (Bayesian FDR selection), [`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md) (pooled permutation null) |
| Matrix and edge utilities | [`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md), [`threshold_to_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/threshold_to_adjacency.md), [`diff_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/diff_adjacency.md), [`upper_tri_to_matrix()`](https://mljaniczek.github.io/ggmtools/reference/upper_tri_to_matrix.md), [`diff_edge_table()`](https://mljaniczek.github.io/ggmtools/reference/diff_edge_table.md) |
| Data | `brca_pi3k`: TCGA BRCA PI3K-pathway expression, Normal vs Tumor (100 x 28 each) |

## Edge tables

Every edge table has the columns `i`, `j`, `node1`, `node2` (one row per
pair i \< j, in [`upper.tri()`](https://rdrr.io/r/base/lower.tri.html)
order), followed by method-specific columns.
[`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md)
adds `group` and `pip`. Differential-edge tables (e.g.
`spikeyglass::diff_edges()`) add `est_g1`, `est_g2`, `diff` and
inference columns, using `lower`/`upper` for intervals.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("mljaniczek/ggmtools")
```
