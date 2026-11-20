# Convert a precision matrix to partial correlations

Computes \\\rho\_{ij} = -\Omega\_{ij} /
\sqrt{\Omega\_{ii}\Omega\_{jj}}\\, with unit diagonal.

## Usage

``` r
precision_to_pcor(Omega)
```

## Arguments

- Omega:

  A square precision (inverse covariance) matrix.

## Value

A symmetric matrix of partial correlations with the dimnames of `Omega`.
If any diagonal entry of `Omega` is not finite and positive, a matrix of
`NA` is returned with a warning.

## Examples

``` r
Omega <- matrix(c(2, -0.8, 0, -0.8, 2, -0.5, 0, -0.5, 1.5), 3, 3)
precision_to_pcor(Omega)
#>      [,1]      [,2]      [,3]
#> [1,]  1.0 0.4000000 0.0000000
#> [2,]  0.4 1.0000000 0.2886751
#> [3,]  0.0 0.2886751 1.0000000
```
