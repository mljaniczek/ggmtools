# Relative Frobenius error of a matrix estimate

Computes \\\\\hat M - M\\\_F / \\M\\\_F\\ over all entries, including
the diagonal.

## Usage

``` r
relative_frobenius(est_mat, truth_mat)
```

## Arguments

- est_mat:

  Estimated matrix (e.g. a precision matrix).

- truth_mat:

  True matrix of the same dimensions.

## Value

A number; `NA` if `truth_mat` is all zero or either matrix has missing
entries.

## See also

[`kl_divergence()`](https://mljaniczek.github.io/ggmtools/reference/kl_divergence.md)

## Examples

``` r
Omega <- diag(3); Omega[1, 2] <- Omega[2, 1] <- 0.4
relative_frobenius(diag(3), Omega)
#> [1] 0.3104602
```
