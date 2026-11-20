# Overlap of the top-k edges of two score vectors

Counts how many edges are among the `k` highest-scoring edges of both
`score1` and `score2`, e.g. two methods' differential scores. Ties are
handled without depending on the order of the edges: when the k-th place
falls inside a block of tied scores, each edge of that block is in the
top k with probability (places left) / (block size), and the result is
the expected overlap under random tie-breaking, \\\sum_e P(e \in A_k)
P(e \in B_k)\\. Without ties at the boundary this is the ordinary
overlap.

## Usage

``` r
topk_overlap(score1, score2, k)
```

## Arguments

- score1, score2:

  Numeric vectors of edge scores of the same length, or p x p matrices
  (the upper triangle is used). Edges with a missing score in either are
  dropped.

- k:

  Integer vector of list sizes.

## Value

A numeric vector of the same length as `k` with the (expected) overlap.
Values of `k` larger than the number of edges are truncated.

## Details

Sorting with [`order()`](https://rdrr.io/r/base/order.html) and taking
the first `k` breaks ties by position. For sparse estimators, which tie
most edges at 0, that inflates the overlap of two methods once `k`
exceeds their number of non-zero edges, because both take the same
low-index edges first.

## Examples

``` r
a <- c(0.9, 0.8, 0, 0, 0, 0)
b <- c(0, 0.7, 0.6, 0, 0, 0)
topk_overlap(a, b, k = 1:6)
#> [1] 0.0000 1.0000 1.6875 2.7500 4.1875 6.0000
```
