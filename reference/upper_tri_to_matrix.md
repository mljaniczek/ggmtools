# Symmetric matrix from an upper triangle

Inverse of `m[upper.tri(m)]`: rebuilds a symmetric p x p matrix from its
upper triangle in [`upper.tri()`](https://rdrr.io/r/base/lower.tri.html)
order (the order of
[`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md)).
Useful for methods that return one value per edge, e.g. a vector of
partial correlations.

## Usage

``` r
upper_tri_to_matrix(x, varnames = NULL, diag = 1)
```

## Arguments

- x:

  Numeric vector of length \\p(p - 1)/2\\.

- varnames:

  Optional character vector of length `p` used as dimnames.

- diag:

  Value for the diagonal. Default 1 (partial correlations).

## Value

A symmetric p x p matrix.

## See also

[`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md)

## Examples

``` r
upper_tri_to_matrix(c(0.5, 0.1, -0.3), c("a", "b", "c"))
#>     a    b    c
#> a 1.0  0.5  0.1
#> b 0.5  1.0 -0.3
#> c 0.1 -0.3  1.0
```
