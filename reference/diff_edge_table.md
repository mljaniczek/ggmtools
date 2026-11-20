# Table of differential edges between two estimates

Lists the edges whose estimates differ by more than `threshold` between
two groups, in the shared edge-table format.

## Usage

``` r
diff_edge_table(mat1, mat2, threshold = 0.05, varnames = NULL)
```

## Arguments

- mat1, mat2:

  Symmetric p x p matrices of edge estimates for groups 1 and 2 (e.g.
  partial correlations).

- threshold:

  Non-negative number; edges with \\\|mat1 - mat2\| \> threshold\\ are
  kept. Default 0.05. Use `-Inf` to keep every edge.

- varnames:

  Optional variable names; defaults to the row names of `mat1`, then
  `V1, ..., Vp`.

## Value

A data frame with columns `i`, `j`, `node1`, `node2` (as in
[`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md)),
`est_g1`, `est_g2` and `diff` (`est_g1 - est_g2`), in
[`upper.tri()`](https://rdrr.io/r/base/lower.tri.html) order. Edges with
a missing difference are dropped. Sort by `-abs(diff)` to list the
largest changes first.

## See also

[`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md),
[`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md)

## Examples

``` r
P1 <- matrix(c(1, 0.4, 0, 0.4, 1, 0.3, 0, 0.3, 1), 3, 3)
P2 <- matrix(c(1, 0.1, 0, 0.1, 1, 0.28, 0, 0.28, 1), 3, 3)
diff_edge_table(P1, P2)
#>   i j node1 node2 est_g1 est_g2 diff
#> 1 1 2    V1    V2    0.4    0.1  0.3
```
