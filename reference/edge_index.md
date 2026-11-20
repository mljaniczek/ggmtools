# Edge index in the shared edge-table format

The first four columns of every edge table in the shared format: one row
per variable pair \\i \< j\\, in
[`upper.tri()`](https://rdrr.io/r/base/lower.tri.html) order
(column-major, i.e. (1, 2), (1, 3), (2, 3), (1, 4), ...).

## Usage

``` r
edge_index(p, varnames = NULL)
```

## Arguments

- p:

  Number of variables.

- varnames:

  Optional character vector of length `p`. Defaults to `V1, ..., Vp`.

## Value

A data frame with columns `i`, `j`, `node1`, `node2`.

## Examples

``` r
edge_index(4, c("a", "b", "c", "d"))
#>   i j node1 node2
#> 1 1 2     a     b
#> 2 1 3     a     c
#> 3 2 3     b     c
#> 4 1 4     a     d
#> 5 2 4     b     d
#> 6 3 4     c     d
```
