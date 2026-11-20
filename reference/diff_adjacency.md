# Structural differential edges between two graphs

An edge is differential when it is present in exactly one of the two
graphs (exclusive or). For differential edges defined by the size of the
change instead, use
`threshold_to_adjacency(theta1 - theta2, threshold)`.

## Usage

``` r
diff_adjacency(adj1, adj2)
```

## Arguments

- adj1, adj2:

  Square matrices of the same dimensions; non-zero entries are edges.

## Value

An integer 0/1 matrix with zero diagonal and the dimnames of `adj1`.
Missing entries stay `NA`.

## See also

[`threshold_to_adjacency()`](https://mljaniczek.github.io/ggmtools/reference/threshold_to_adjacency.md)

## Examples

``` r
a1 <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3)
a2 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
diff_adjacency(a1, a2)
#>      [,1] [,2] [,3]
#> [1,]    0    0    1
#> [2,]    0    0    1
#> [3,]    1    1    0
```
