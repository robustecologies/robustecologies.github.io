# Print robust Bayesian sensitivity analysis results

Concise one-screen summary of a `robust_bayes` object showing key
results from each module.

## Usage

``` r
# S3 method for class 'robust_bayes'
print(x, ...)
```

## Arguments

- x:

  An object of class `robust_bayes`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a RobustBayes
prior-sensitivity analysis produced by
[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Berger, J. O., & Berliner, L. M. (1986). Robust Bayes and empirical
Bayes analysis with epsilon-contaminated priors. *Annals of Statistics*,
14(2), 461-486.
[doi:10.1214/aos/1176349933](https://doi.org/10.1214/aos/1176349933)

## See also

[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md),
[`summary.robust_bayes`](https://robustecologies.github.io/lucifer/reference/summary.robust_bayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.robust_bayes
} # }
```
