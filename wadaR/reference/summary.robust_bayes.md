# Summarize robust Bayesian sensitivity analysis results

Extended output with full tables for each module including complete
sensitivity scores, divergence matrices, conflict classifications,
influence statistics, and weight comparisons.

## Usage

``` r
# S3 method for class 'robust_bayes'
summary(object, ...)
```

## Arguments

- object:

  An object of class `robust_bayes`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns `object`.

## Details

Produces a tabular summary of a RobustBayes prior-sensitivity analysis
produced by
[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Berger, J. O., & Berliner, L. M. (1986). Robust Bayes and empirical
Bayes analysis with epsilon-contaminated priors. *Annals of Statistics*,
14(2), 461-486.
[doi:10.1214/aos/1176349933](https://doi.org/10.1214/aos/1176349933)

## See also

[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md),
[`print.robust_bayes`](https://robustecologies.github.io/lucifer/reference/print.robust_bayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.robust_bayes
} # }
```
