# Summarize data cloning MLE results

Extended output with full MLE table, per-K diagnostic tables, eigenvalue
ratios, MSE/R-squared, scaled variance convergence, and ANOVA
estimability results (when multiple priors are used).

## Usage

``` r
# S3 method for class 'data_cloning'
summary(object, ...)
```

## Arguments

- object:

  An object of class `data_cloning`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns `object`.

## Details

Produces a tabular summary of a data-cloning identifiability analysis
produced by
[`DataCloning`](https://robustecologies.github.io/lucifer/reference/DataCloning.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Lele, S. R., Dennis, B., & Lutscher, F. (2007). Data cloning: easy
maximum likelihood estimation for complex ecological models using
Bayesian Markov chain Monte Carlo methods. *Ecology Letters*, 10(7),
551-563.
[doi:10.1111/j.1461-0248.2007.01047.x](https://doi.org/10.1111/j.1461-0248.2007.01047.x)

## See also

[`DataCloning`](https://robustecologies.github.io/lucifer/reference/DataCloning.md),
[`print.data_cloning`](https://robustecologies.github.io/lucifer/reference/print.data_cloning.md),
[`plot.data_cloning`](https://robustecologies.github.io/lucifer/reference/plot.data_cloning.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.data_cloning
} # }
```
