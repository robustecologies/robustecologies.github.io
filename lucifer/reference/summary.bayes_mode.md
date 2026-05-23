# Summary method for Bayesian mode inference

Provides an extended summary including MCMC diagnostics, effective
number of components, mixture fit quality, and detailed mode location
posterior distributions.

## Usage

``` r
# S3 method for class 'bayes_mode'
summary(object, digits = 4, ...)
```

## Arguments

- object:

  An object of class `bayes_mode` or `bayes_mode_multi`.

- digits:

  Number of decimal digits. Default 4.

- ...:

  Currently unused.

## Value

Invisibly returns a list of class `summary.bayes_mode`.

## Details

Produces a tabular summary of a Bayesian mode inference fit produced by
[`BayesMode`](https://robustecologies.github.io/lucifer/reference/BayesMode.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Cross, J. L., Hoogerheide, L., Ulker, P., & van Dijk, H. K. (2024).
Sparse finite mixtures for modal inference. *Economics Letters*, 235,
111551.
[doi:10.1016/j.econlet.2024.111551](https://doi.org/10.1016/j.econlet.2024.111551)

## See also

[`plot.bayes_mode`](https://robustecologies.github.io/lucifer/reference/plot.bayes_mode.md),
[`print.bayes_mode`](https://robustecologies.github.io/lucifer/reference/print.bayes_mode.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.bayes_mode
} # }
```
