# Partial correlations from a joint graphical lasso fit

[`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md)
method for fits from
[`JGL::JGL()`](https://rdrr.io/pkg/JGL/man/JGL.html) (class `"jgl"`;
Danaher, Wang & Witten 2014): converts each group's estimated precision
matrix with
[`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md).
It only reads `fit$theta`, so ggmtools does not depend on JGL. (JGL was
archived on CRAN in July 2025; it can be installed from the CRAN archive
or with `remotes::install_github("cran/JGL")`.)

## Usage

``` r
# S3 method for class 'jgl'
pcor_estimate(fit, varnames = NULL, groupnames = NULL, ...)
```

## Arguments

- fit:

  A `"jgl"` object from
  [`JGL::JGL()`](https://rdrr.io/pkg/JGL/man/JGL.html), fitted with
  `return.whole.theta = TRUE`. With the default
  `return.whole.theta = FALSE`, JGL drops variables that are not
  connected to any other, and the method stops with an error.

- varnames:

  Optional variable names. Defaults to the dimnames of the estimates
  (JGL keeps the column names of the data), then `V1, ..., Vp`.

- groupnames:

  Optional group names. JGL does not store them, so they default to
  `Group_1, ..., Group_K`; pass `names(data_list)` to use the names of
  the data list.

- ...:

  Ignored.

## Value

A named list of K symmetric p x p partial-correlation matrices, as
described in
[`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md).
These are penalized point estimates: entries the penalty sets to zero
are exactly 0.

## References

Danaher, P., Wang, P. & Witten, D.M. (2014). The joint graphical lasso
for inverse covariance estimation across multiple classes. *Journal of
the Royal Statistical Society: Series B*, 76(2), 373-397.

## See also

[`ggm_bic()`](https://mljaniczek.github.io/ggmtools/reference/ggm_bic.md)
for tuning the penalties.

## Examples

``` r
# A minimal object with the structure JGL::JGL() returns
theta <- list(diag(3), matrix(c(2, -0.5, 0, -0.5, 2, 0, 0, 0, 1), 3, 3))
fit <- structure(list(theta = theta, connected = rep(TRUE, 3)), class = "jgl")
pcor_estimate(fit, varnames = c("a", "b", "c"), groupnames = c("A", "B"))
#> $A
#>   a b c
#> a 1 0 0
#> b 0 1 0
#> c 0 0 1
#> 
#> $B
#>      a    b c
#> a 1.00 0.25 0
#> b 0.25 1.00 0
#> c 0.00 0.00 1
#> 
```
