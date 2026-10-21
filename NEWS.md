# ggmtools (development version)

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
