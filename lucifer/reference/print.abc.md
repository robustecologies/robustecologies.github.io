# Print method for ABC objects

Prints a concise summary of an ABC fit.

## Usage

``` r
# S3 method for class 'abc'
print(x, ...)
```

## Arguments

- x:

  An object of class `abc`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of an approximate Bayesian
computation fit produced by
[`ABC`](https://robustecologies.github.io/lucifer/reference/ABC.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Beaumont, M. A., Zhang, W., & Balding, D. J. (2002). Approximate
Bayesian computation in population genetics. *Genetics*, 162(4),
2025-2035.
[doi:10.1093/genetics/162.4.2025](https://doi.org/10.1093/genetics/162.4.2025)

## See also

[`plot.abc`](https://robustecologies.github.io/lucifer/reference/plot.abc.md),
[`summary.abc`](https://robustecologies.github.io/lucifer/reference/summary.abc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.abc
} # }
```
