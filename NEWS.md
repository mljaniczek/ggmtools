# ggmtools 0.0.0.9000

## Evaluation and inference helpers from the DiNA simulation

Method-agnostic helpers moved from the DiNA simulation
(`dina_simulation/04_calc_metrics.R`, `00_helper_functions.R`,
`tcga_demo/utils.R`), with the overlapping versions in spikeyglass and
multiGGMr merged into one.

* New: `pr_auc()` (PR curve and AUPR; `method = "trapezoid"`,
  `"average_precision"` or `"interpolated"`), `topk_overlap()`,
  `bayes_fdr()`, `pooled_perm_pvalues()`, `relative_frobenius()`,
  `kl_divergence()`, `ggm_bic()`, `threshold_to_adjacency()`,
  `diff_adjacency()`, `upper_tri_to_matrix()` and `diff_edge_table()`.
* `confusion_at_threshold()` also returns `FDR` and `MCC`.
* `confusion_at_threshold()`, `roc_auc()` and `pr_auc()` accept vectors of
  edge scores and truths as well as p x p matrices.
* `stats` added to Imports.

Old to new names:

| Old | New |
|---|---|
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
| DiNA `bayes_fdr_call()`, `.bayes_fdr_call()`; multiGGMr `bayes_fdr()` | `bayes_fdr(prob, alpha)` |
| DiNA `pooled_perm_pvalues()`; spikeyglass `.pooled_perm_pvalues()` | `pooled_perm_pvalues(observed, null)` |

Behaviour differences from the old versions:

* `bayes_fdr()` keeps tied probabilities together: the cutoff is placed only
  between distinct values, so the selection no longer depends on the input
  order. If a whole tie block would push the expected FDR above `alpha`,
  none of it is selected. Example: `prob = c(1, 1, 0.9, 0.9, 0.9)`,
  `alpha = 0.05` used to select 4 edges (two of the three tied 0.9s, chosen
  by position); it now selects 2. Results without ties at the cutoff are
  unchanged.
* `pooled_perm_pvalues()` is spikeyglass's version: a null value counts as
  an exceedance when it is within a relative tolerance
  (`sqrt(.Machine$double.eps) * max(1, |x|)`) of the observed score. DiNA
  subtracted an absolute `.Machine$double.eps`, which vanishes for scores
  `>= 2`, so exact ties were then not counted. Example: observed 3 against
  the null `c(1, 2, 3, 3, 3, 4)` gave p = 2/7 and now gives 5/7. Without
  ties the p-values are unchanged. It returns `pvalue` and `p_adj` (DiNA:
  `pvalues`, `p_adj`, `upper_tri_idx`; use `edge_index(p)` for the indices).
* `kl_divergence()` returns the Kullback-Leibler divergence
  KL(N(0, Sigma) || N(0, theta_hat^-1)) =
  (tr(Sigma theta_hat) - log|Sigma theta_hat| - p) / 2. DiNA's
  `compute_kl_divergence()` omitted the 1/2 (it is Stein's loss), so its
  values are twice `kl_divergence()`; its comment also had the direction
  reversed. Rankings are unchanged. A singular `theta_hat` gives `NA`
  instead of `Inf`.
* `ggm_bic()` equals `vBIC()` for valid input; it returns `NA` when an
  estimated precision matrix has missing values or a zero or negative
  determinant (`vBIC()` used `log|det|` and ignored the sign). It takes the
  list of precision matrices, not the fitted JGL object.
* `pr_auc()` with the default `method = "trapezoid"` equals DiNA's
  tie-corrected `compute_aupr()`. `"average_precision"` and
  `"interpolated"` (Davis-Goadrich) are new; the trapezoid is optimistic
  when a large tied block (e.g. edges at exactly 0) has lower precision
  than the edges above it.
* `precision_to_pcor()` symmetrizes the result and returns `NA` with a
  warning for a non-positive diagonal; DiNA's `pcor_from_precision()` kept
  asymmetry and returned `NaN`.
* `diff_edge_table()` returns the shared edge-table columns `i`, `j`,
  `node1`, `node2`, `est_g1`, `est_g2`, `diff` in `upper.tri()` order. DiNA's
  `extract_diff_edges()` returned `gene1`, `gene2`, `pcor_group1`,
  `pcor_group2`, `diff`, sorted by `-abs(diff)`.
* `topk_overlap()` counts the overlap of two top-k lists with ties broken
  at random in expectation. Taking `order(...)[1:k]`, as DiNA's
  `plot_topk_overlap()` does, breaks ties by position, which inflates the
  overlap of sparse methods once `k` exceeds their non-zero edges.

## Initial version

* Initial version. Collects the functions previously duplicated in
  spikeyglass and multiGGMr (`precision_to_pcor()`, `confusion_at_threshold()`,
  `roc_auc()`, `plot_roc()`), with one signature, plus the shared generics
  `pip_table()` and `pcor_estimate()`, the edge-table helpers `edge_index()`
  and `as_pip_table()`, and the `brca_pi3k` dataset.
* `roc_auc()` handles tied scores correctly: the curve steps only at
  distinct scores, so the AUC equals the Mann-Whitney statistic. (multiGGMr's
  previous version accumulated one edge at a time, so tied edges were ordered
  by matrix position.) It returns `auc` (lower case).
* `confusion_at_threshold()` returns a list including `precision` and `F1`.
  Rates with a zero denominator are `NA`.
