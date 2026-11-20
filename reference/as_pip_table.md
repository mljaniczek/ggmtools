# Build a long PIP table from per-group PIP matrices

Helper for implementing
[`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md)
methods.

## Usage

``` r
as_pip_table(pip_list, varnames = NULL, groupnames = NULL)
```

## Arguments

- pip_list:

  A list of K symmetric p x p PIP matrices (one per group; repeat the
  same matrix if the PIP is shared across groups).

- varnames:

  Optional variable names; defaults to the dimnames of the first matrix,
  then `V1, ..., Vp`.

- groupnames:

  Optional group names; defaults to `names(pip_list)`, then
  `Group_1, ..., Group_K`.

## Value

A data frame in the format described in
[`pip_table()`](https://mljaniczek.github.io/ggmtools/reference/pip_table.md).

## Examples

``` r
m <- matrix(c(0, .9, .1, .9, 0, .4, .1, .4, 0), 3, 3)
as_pip_table(list(A = m, B = m))
#>   i j node1 node2 group pip
#> 1 1 2    V1    V2     A 0.9
#> 2 1 3    V1    V3     A 0.1
#> 3 2 3    V2    V3     A 0.4
#> 4 1 2    V1    V2     B 0.9
#> 5 1 3    V1    V3     B 0.1
#> 6 2 3    V2    V3     B 0.4
```
