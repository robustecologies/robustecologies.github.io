# Print data cloning MLE results

Concise one-screen summary of a `data_cloning` object showing the MLE,
standard errors, Wald confidence intervals, and overall convergence
status.

## Usage

``` r
# S3 method for class 'data_cloning'
print(x, ...)
```

## Arguments

- x:

  An object of class `data_cloning`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a data-cloning
identifiability analysis produced by
[`DataCloning`](https://robustecologies.github.io/lucifer/reference/DataCloning.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Lele, S. R., Dennis, B., & Lutscher, F. (2007). Data cloning: easy
maximum likelihood estimation for complex ecological models using
Bayesian Markov chain Monte Carlo methods. *Ecology Letters*, 10(7),
551-563.
[doi:10.1111/j.1461-0248.2007.01047.x](https://doi.org/10.1111/j.1461-0248.2007.01047.x)

## See also

[`DataCloning`](https://robustecologies.github.io/lucifer/reference/DataCloning.md),
[`summary.data_cloning`](https://robustecologies.github.io/lucifer/reference/summary.data_cloning.md),
[`plot.data_cloning`](https://robustecologies.github.io/lucifer/reference/plot.data_cloning.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.data_cloning
} # }
```
