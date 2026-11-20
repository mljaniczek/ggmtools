# Threshold a matrix to an adjacency matrix

Declares an edge where \\\|m\_{ij}\| \> threshold\\, e.g. to turn an
estimated precision or partial-correlation matrix into a graph. Applied
to a difference `theta1 - theta2`, it gives magnitude-based differential
edges.

## Usage

``` r
threshold_to_adjacency(mat, threshold = 0.01)
```

## Arguments

- mat:

  A square numeric matrix.

- threshold:

  Non-negative number. Default 0.01.

## Value

An integer 0/1 matrix with zero diagonal and the dimnames of `mat`.
Missing entries stay `NA`.

## See also

[`diff_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/diff_adjacency.md),
[`confusion_at_threshold()`](https://mljaniczek.github.io/ggmtools/reference/confusion_at_threshold.md)

## Examples

``` r
Omega <- matrix(c(2, -0.8, 0.005, -0.8, 2, -0.5, 0.005, -0.5, 1.5), 3, 3)
threshold_to_adjacency(Omega)
#>      [,1] [,2] [,3]
#> [1,]    0    1    0
#> [2,]    1    0    1
#> [3,]    0    1    0
```
