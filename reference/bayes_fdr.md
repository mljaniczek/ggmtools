# Bayesian false discovery rate selection

Selects the largest set of hypotheses whose posterior expected false
discovery rate is at most `alpha`, using the direct posterior
probability approach (Newton et al. 2004). Hypotheses are ranked by
their posterior probability `prob`; the expected FDR of the top \\k\\ is
the mean of \\1 - prob\\ over them.

## Usage

``` r
bayes_fdr(prob, alpha = 0.05)
```

## Arguments

- prob:

  Numeric vector of posterior probabilities (e.g. posterior inclusion
  probabilities), or a symmetric matrix of edge probabilities (only the
  upper triangle is used).

- alpha:

  Target posterior expected FDR. Default 0.05.

## Value

A logical vector, or a symmetric logical matrix with `FALSE` diagonal,
marking the selected hypotheses. Attributes: `expected_fdr` (posterior
expected FDR of the selected set, 0 if none), `n_selected`, `cutoff`
(smallest selected probability, or `NA`) and `alpha`.

## Details

The guarantee is only as good as the calibration of `prob`. For
differential edges, pass a probability of a real difference such as
\\P(\|D\| \> c)\\ or \\P(g_k \ne g_m)\\, not a sign probability.

Tied probabilities are selected or rejected together: the cutoff is only
placed between distinct values of `prob`, so the selection does not
depend on the order of the input. If including a whole block of ties
would push the expected FDR above `alpha`, none of the block is
selected.

Missing probabilities are never selected and do not count towards the
expected FDR.

## References

Newton, M.A., Noueiry, A., Sarkar, D. & Ahlquist, P. (2004). Detecting
differential gene expression with a semiparametric hierarchical mixture
method. *Biostatistics*, 5(2), 155-176.

## See also

[`pooled_perm_pvalues()`](https://mljaniczek.github.io/ggmtools/reference/pooled_perm_pvalues.md)
for permutation-based selection.

## Examples

``` r
bayes_fdr(c(0.99, 0.97, 0.9, 0.6, 0.2), alpha = 0.05)
#> [1]  TRUE  TRUE  TRUE FALSE FALSE
#> attr(,"expected_fdr")
#> [1] 0.04666667
#> attr(,"n_selected")
#> [1] 3
#> attr(,"cutoff")
#> [1] 0.9
#> attr(,"alpha")
#> [1] 0.05

# Tied probabilities are kept together
bayes_fdr(c(1, 1, 0.9, 0.9, 0.9), alpha = 0.05)
#> [1]  TRUE  TRUE FALSE FALSE FALSE
#> attr(,"expected_fdr")
#> [1] 0
#> attr(,"n_selected")
#> [1] 2
#> attr(,"cutoff")
#> [1] 1
#> attr(,"alpha")
#> [1] 0.05
```
