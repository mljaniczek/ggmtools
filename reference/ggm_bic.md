# BIC for multiple Gaussian graphical models

A BIC-type criterion for tuning joint estimators of K precision matrices
(e.g. the penalties of the joint graphical lasso), summed over groups:
\$\$\sum_k \left\[\mathrm{tr}(S_k\hat\Theta_k) - \log\|\hat\Theta_k\| +
\frac{\log n_k}{n_k} E_k\right\],\$\$ where \\S_k\\ is the sample
covariance of group \\k\\ (denominator \\n_k - 1\\), \\n_k\\ its sample
size and \\E_k\\ the number of entries of \\\hat\Theta_k\\ with
\\\|\hat\theta\_{ij}\| \ge\\ `threshold`. \\E_k\\ counts every such
entry of the matrix, including the diagonal and both triangles. This is
the per-observation form used in the Koyejo lab joint graphical lasso
code (`vBIC()`); smaller is better.

## Usage

``` r
ggm_bic(data_list, theta_list, threshold = 0.001)
```

## Arguments

- data_list:

  List of K data matrices (\\n_k\\ x p).

- theta_list:

  List of K estimated p x p precision matrices.

- threshold:

  Entries with absolute value `>= threshold` count as non-zero. Default
  0.001.

## Value

A number; `NA` if any \\\hat\Theta_k\\ has missing values or a
non-positive determinant, or the data have missing values.

## References

Danaher, P., Wang, P. & Witten, D.M. (2014). The joint graphical lasso
for inverse covariance estimation across multiple classes. *Journal of
the Royal Statistical Society: Series B*, 76(2), 373-397.

## Examples

``` r
set.seed(1)
X <- list(matrix(rnorm(100), 25, 4), matrix(rnorm(80), 20, 4))
theta <- lapply(X, function(x) solve(cov(x)))
ggm_bic(X, theta)
#> [1] 10.12662
```
