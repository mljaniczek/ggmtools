# Package index

## Edge tables and generics

The shared long edge-table format, and generics that packages implement
for their fitted models.

- [`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md)
  : Edge index in the shared edge-table format
- [`as_pip_table()`](https://mljaniczek.github.io/ggmtools/reference/as_pip_table.md)
  : Build a long PIP table from per-group PIP matrices
- [`diff_edge_table()`](https://mljaniczek.github.io/ggmtools/reference/diff_edge_table.md)
  : Table of differential edges between two estimates
- [`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md)
  : Posterior inclusion probabilities as a long edge table
- [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md)
  : Point estimates of partial correlations
- [`pcor_estimate(`*`<jgl>`*`)`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.jgl.md)
  : Partial correlations from a joint graphical lasso fit

## Evaluation against a known graph

- [`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md)
  : ROC curve and AUC for edge recovery
- [`pr_auc()`](https://mljaniczek.github.io/ggmtools/reference/pr_auc.md)
  : Precision-recall curve and AUPR for edge recovery
- [`plot_roc()`](https://mljaniczek.github.io/ggmtools/reference/plot_roc.md)
  : Plot an ROC curve
- [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md)
  : Confusion metrics for edge selection at a threshold
- [`topk_overlap()`](https://mljaniczek.github.io/ggmtools/reference/topk_overlap.md)
  : Overlap of the top-k edges of two score vectors

## Estimation error and tuning

- [`relative_frobenius()`](https://mljaniczek.github.io/ggmtools/reference/relative_frobenius.md)
  : Relative Frobenius error of a matrix estimate
- [`kl_divergence()`](https://mljaniczek.github.io/ggmtools/reference/kl_divergence.md)
  : Kullback-Leibler divergence of an estimated Gaussian graphical model
- [`ggm_bic()`](https://mljaniczek.github.io/ggmtools/reference/ggm_bic.md)
  : BIC for multiple Gaussian graphical models

## Edge selection

- [`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md)
  : Bayesian false discovery rate selection
- [`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md)
  : Empirical p-values against a pooled permutation null

## Matrix utilities

- [`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md)
  : Convert a precision matrix to partial correlations
- [`threshold_to_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/threshold_to_adjacency.md)
  : Threshold a matrix to an adjacency matrix
- [`diff_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/diff_adjacency.md)
  : Structural differential edges between two graphs
- [`upper_tri_to_matrix()`](https://mljaniczek.github.io/ggmtools/reference/upper_tri_to_matrix.md)
  : Symmetric matrix from an upper triangle

## Data

- [`brca_pi3k`](https://mljaniczek.github.io/ggmtools/reference/brca_pi3k.md)
  : TCGA breast cancer PI3K pathway expression, Normal vs Tumor
