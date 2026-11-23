# Changelog

## ggmtools 0.0.0.9000

### Evaluation and inference helpers from the DiNA simulation

Method-agnostic helpers moved from the DiNA simulation
(`dina_simulation/04_calc_metrics.R`, `00_helper_functions.R`,
`tcga_demo/utils.R`), with the overlapping versions in spikeyglass and
multiGGMr merged into one.

- New:
  [`pr_auc()`](https://mljaniczek.github.io/ggmtools/reference/pr_auc.md)
  (PR curve and AUPR; `method = "trapezoid"`, `"average_precision"` or
  `"interpolated"`),
  [`topk_overlap()`](https://mljaniczek.github.io/ggmtools/reference/topk_overlap.md),
  [`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md),
  [`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md),
  [`relative_frobenius()`](https://mljaniczek.github.io/ggmtools/reference/relative_frobenius.md),
  [`kl_divergence()`](https://mljaniczek.github.io/ggmtools/reference/kl_divergence.md),
  [`ggm_bic()`](https://mljaniczek.github.io/ggmtools/reference/ggm_bic.md),
  [`threshold_to_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/threshold_to_adjacency.md),
  [`diff_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/diff_adjacency.md),
  [`upper_tri_to_matrix()`](https://mljaniczek.github.io/ggmtools/reference/upper_tri_to_matrix.md)
  and
  [`diff_edge_table()`](https://mljaniczek.github.io/ggmtools/reference/diff_edge_table.md).
- [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md)
  also returns `FDR` and `MCC`.
- [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md),
  [`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md)
  and
  [`pr_auc()`](https://mljaniczek.github.io/ggmtools/reference/pr_auc.md)
  accept vectors of edge scores and truths as well as p x p matrices.
- `stats` added to Imports.
- New
  [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md)
  method for joint graphical lasso fits from
  [`JGL::JGL()`](https://rdrr.io/pkg/JGL/man/JGL.html) (class `"jgl"`).
  It only reads `fit$theta`, so ggmtools does not depend on JGL. Fits
  must use `return.whole.theta = TRUE`.

Old to new names:

| Old | New |
|----|----|
| DiNA `compute_auc(truth_vec, score_vec)` | `roc_auc(score_vec, truth_vec)$auc` |
| DiNA `compute_aupr(truth_vec, score_vec)` | `pr_auc(score_vec, truth_vec)$auc` |
| DiNA `compute_tpr_fpr(truth_adj, est_adj)` | `confusion_at_threshold(est_adj, truth_adj)` |
| DiNA `compute_fdr(TP, FP)`, `compute_mcc(TP, FP, TN, FN)` | `confusion_at_threshold(...)$FDR`, `$MCC` |
| DiNA `compute_relative_frobenius(est, truth)` | `relative_frobenius(est, truth)` |
| DiNA `compute_kl_divergence(theta_hat, sigma_true)` | `2 * kl_divergence(theta_hat, sigma_true)` |
| DiNA `threshold_to_adjacency(theta, threshold)` | `threshold_to_adjacency(theta, threshold)` |
| DiNA XOR differential truth (`.xor_truth()` builds it) | `diff_adjacency(adj1, adj2)` |
| DiNA `pcor_from_precision(Theta)` | `precision_to_pcor(Theta)` |
| DiNA `vBIC(X, est_graph, thr)` | `ggm_bic(X, est_graph$theta, threshold = thr)` |
| DiNA `vec_to_pcor_mat(vec, p, gene_names)` | `upper_tri_to_matrix(vec, gene_names)` |
| DiNA `extract_diff_edges(pcor1, pcor2, threshold, gene_names)` | `diff_edge_table(pcor1, pcor2, threshold, gene_names)` |
| DiNA `bayes_fdr_call()`, `.bayes_fdr_call()`; multiGGMr [`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md) | `bayes_fdr(prob, alpha)` |
| DiNA [`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md); spikeyglass `.pooled_perm_pvalues()` | `pooled_perm_pvalues(observed, null)` |

Behaviour differences from the old versions:

- [`bayes_fdr()`](https://mljaniczek.github.io/ggmtools/reference/bayes_fdr.md)
  keeps tied probabilities together: the cutoff is placed only between
  distinct values, so the selection no longer depends on the input
  order. If a whole tie block would push the expected FDR above `alpha`,
  none of it is selected. Example: `prob = c(1, 1, 0.9, 0.9, 0.9)`,
  `alpha = 0.05` used to select 4 edges (two of the three tied 0.9s,
  chosen by position); it now selects 2. Results without ties at the
  cutoff are unchanged.
- [`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md)
  is spikeyglass’s version: a null value counts as an exceedance when it
  is within a relative tolerance
  (`sqrt(.Machine$double.eps) * max(1, |x|)`) of the observed score.
  DiNA subtracted an absolute `.Machine$double.eps`, which vanishes for
  scores `>= 2`, so exact ties were then not counted. Example: observed
  3 against the null `c(1, 2, 3, 3, 3, 4)` gave p = 2/7 and now gives
  5/7. Without ties the p-values are unchanged. It returns `pvalue` and
  `p_adj` (DiNA: `pvalues`, `p_adj`, `upper_tri_idx`; use
  `edge_index(p)` for the indices).
- [`kl_divergence()`](https://mljaniczek.github.io/ggmtools/reference/kl_divergence.md)
  returns the Kullback-Leibler divergence KL(N(0, Sigma) \|\| N(0,
  theta_hat^-1)) = (tr(Sigma theta_hat) - log\|Sigma theta_hat\| - p)
  / 2. DiNA’s `compute_kl_divergence()` omitted the 1/2 (it is Stein’s
  loss), so its values are twice
  [`kl_divergence()`](https://mljaniczek.github.io/ggmtools/reference/kl_divergence.md);
  its comment also had the direction reversed. Rankings are unchanged. A
  singular `theta_hat` gives `NA` instead of `Inf`.
- [`ggm_bic()`](https://mljaniczek.github.io/ggmtools/reference/ggm_bic.md)
  equals `vBIC()` for valid input; it returns `NA` when an estimated
  precision matrix has missing values or a zero or negative determinant
  (`vBIC()` used `log|det|` and ignored the sign). It takes the list of
  precision matrices, not the fitted JGL object.
- [`pr_auc()`](https://mljaniczek.github.io/ggmtools/reference/pr_auc.md)
  with the default `method = "trapezoid"` equals DiNA’s tie-corrected
  `compute_aupr()`. `"average_precision"` and `"interpolated"`
  (Davis-Goadrich) are new; the trapezoid is optimistic when a large
  tied block (e.g. edges at exactly 0) has lower precision than the
  edges above it.
- [`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md)
  symmetrizes the result and returns `NA` with a warning for a
  non-positive diagonal; DiNA’s `pcor_from_precision()` kept asymmetry
  and returned `NaN`.
- [`diff_edge_table()`](https://mljaniczek.github.io/ggmtools/reference/diff_edge_table.md)
  returns the shared edge-table columns `i`, `j`, `node1`, `node2`,
  `est_g1`, `est_g2`, `diff` in
  [`upper.tri()`](https://rdrr.io/r/base/lower.tri.html) order. DiNA’s
  `extract_diff_edges()` returned `gene1`, `gene2`, `pcor_group1`,
  `pcor_group2`, `diff`, sorted by `-abs(diff)`.
- [`topk_overlap()`](https://mljaniczek.github.io/ggmtools/reference/topk_overlap.md)
  counts the overlap of two top-k lists with ties broken at random in
  expectation. Taking `order(...)[1:k]`, as DiNA’s `plot_topk_overlap()`
  does, breaks ties by position, which inflates the overlap of sparse
  methods once `k` exceeds their non-zero edges.

### Initial version

- Initial version. Collects the functions previously duplicated in
  spikeyglass and multiGGMr
  ([`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md),
  [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md),
  [`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md),
  [`plot_roc()`](https://mljaniczek.github.io/ggmtools/reference/plot_roc.md)),
  with one signature, plus the shared generics
  [`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md)
  and
  [`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md),
  the edge-table helpers
  [`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md)
  and
  [`as_pip_table()`](https://mljaniczek.github.io/ggmtools/reference/as_pip_table.md),
  and the `brca_pi3k` dataset.
- [`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md)
  handles tied scores correctly: the curve steps only at distinct
  scores, so the AUC equals the Mann-Whitney statistic. (multiGGMr’s
  previous version accumulated one edge at a time, so tied edges were
  ordered by matrix position.) It returns `auc` (lower case).
- [`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md)
  returns a list including `precision` and `F1`. Rates with a zero
  denominator are `NA`.
