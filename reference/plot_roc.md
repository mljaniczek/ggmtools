# Plot an ROC curve

Plot an ROC curve

## Usage

``` r
plot_roc(roc_obj, main = "ROC curve", add = FALSE, ...)
```

## Arguments

- roc_obj:

  A list from
  [`roc_auc()`](https://mljaniczek.github.io/ggmtools/reference/roc_auc.md).
  (Older lists with an `AUC` element instead of `auc` are accepted.)

- main:

  Plot title.

- add:

  Logical. If `TRUE`, add the curve to the current plot.

- ...:

  Passed to [`plot`](https://rdrr.io/r/graphics/plot.default.html) or
  [`lines`](https://rdrr.io/r/graphics/lines.html).

## Value

`roc_obj`, invisibly.

## Examples

``` r
truth <- matrix(rbinom(100, 1, 0.3), 10, 10); truth <- truth + t(truth)
plot_roc(roc_auc(matrix(runif(100), 10, 10), truth))
```
