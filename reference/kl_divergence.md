# Kullback-Leibler divergence of an estimated Gaussian graphical model

Computes the Kullback-Leibler divergence from the true zero-mean
Gaussian \\N(0, \Sigma)\\ to the estimated one \\N(0,
\hat\Theta^{-1})\\, \$\$KL =
\frac{1}{2}\left\[\mathrm{tr}(\Sigma\hat\Theta) -
\log\|\Sigma\hat\Theta\| - p\right\].\$\$ Twice this quantity is Stein's
(entropy) loss of \\\hat\Theta\\.

## Usage

``` r
kl_divergence(theta_hat, sigma_true)
```

## Arguments

- theta_hat:

  Estimated p x p precision matrix.

- sigma_true:

  True p x p covariance matrix.

## Value

A non-negative number, or `NA` if \\\Sigma\hat\Theta\\ has a zero or
negative determinant (e.g. `theta_hat` is not positive definite), the
inputs have missing values, or the computation fails.

## See also

[`relative_frobenius()`](https://mljaniczek.github.io/ggmtools/reference/relative_frobenius.md)

## Examples

``` r
Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
kl_divergence(solve(Sigma), Sigma)   # 0
#> [1] 0
kl_divergence(diag(2), Sigma)
#> [1] 0.143841
```
