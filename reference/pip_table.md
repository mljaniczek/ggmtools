# Posterior inclusion probabilities as a long edge table

Generic for extracting edge-level posterior inclusion probabilities
(PIPs) from a fitted multiple-GGM model, in a format shared by all
packages that implement it, so that one downstream script works on any
of their fits.

## Usage

``` r
pip_table(fit, ...)
```

## Arguments

- fit:

  A fitted model object.

- ...:

  Passed to methods (e.g. which ladder step or chain to use).

## Value

A data frame with one row per edge \\i \< j\\ and group, ordered by
group and then in [`upper.tri()`](https://rdrr.io/r/base/lower.tri.html)
order, with columns

- i, j:

  Variable indices, \\i \< j\\.

- node1, node2:

  Variable names.

- group:

  Group name (character).

- pip:

  Posterior probability that the edge is in that group's graph.

Methods document how their PIPs are defined. For example, in spikeyglass
the edge indicator is shared by all groups, so `pip` is identical across
groups; in multiGGMr each group has its own graph.

## See also

[`pcor_estimate()`](https://mljaniczek.github.io/ggmtools/reference/pcor_estimate.md),
[`edge_index()`](https://mljaniczek.github.io/ggmtools/reference/edge_index.md)
