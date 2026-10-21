# ggmtools

Shared tools for R packages that estimate multiple Gaussian graphical
models, currently [spikeyglass](https://github.com/mljaniczek/spikeyglass)
(Bayesian spike-and-slab joint graphical lasso, EM) and
[multiGGMr](https://github.com/mljaniczek/multiGGMr) (Bayesian multiple GGMs
with MRF priors, MCMC).

Both packages import ggmtools and re-export what users call. So the
functions below are the same objects whichever package is attached, and
one downstream script works on fits from either.

| What | Functions |
|---|---|
| Shared edge-table format | `edge_index()`, `as_pip_table()` |
| Generics with methods in each package | `pip_table()` (long table of PIPs), `pcor_estimate()` (one partial-correlation matrix per group) |
| Evaluation | `roc_auc()` (tie-correct ROC/AUC), `plot_roc()`, `confusion_at_threshold()` |
| Utilities | `precision_to_pcor()` |
| Data | `brca_pi3k`: TCGA BRCA PI3K-pathway expression, Normal vs Tumor (100 x 28 each) |

## Edge tables

Every edge table has the columns `i`, `j`, `node1`, `node2` (one row per
pair i < j, in `upper.tri()` order), followed by method-specific columns.
`pip_table()` adds `group` and `pip`. Differential-edge tables (e.g.
`spikeyglass::diff_edges()`) add `est_g1`, `est_g2`, `diff` and inference
columns, using `lower`/`upper` for intervals.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mljaniczek/ggmtools")
```
