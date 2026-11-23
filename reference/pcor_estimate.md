# Point estimates of partial correlations

Generic for extracting one partial-correlation matrix per group from a
fitted multiple-GGM model.

## Usage

``` r
pcor_estimate(fit, ...)
```

## Arguments

- fit:

  A fitted model object.

- ...:

  Passed to methods.

## Value

A named list of K symmetric p x p matrices with unit diagonal and
variable names as dimnames. The names are the group names, or
`"Group_1"`, `"Group_2"`, ... if the groups are unnamed. Methods
document which estimate they return. For example, spikeyglass returns
the MAP estimate and multiGGM the posterior mean.

## See also

[`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md),
[`precision_to_pcor()`](https://mljaniczek.github.io/ggmtools/reference/precision_to_pcor.md)
