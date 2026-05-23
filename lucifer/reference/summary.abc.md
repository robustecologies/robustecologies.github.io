# Summary method for ABC objects

Provides an extended summary of an ABC fit.

## Usage

``` r
# S3 method for class 'abc'
summary(object, ...)
```

## Arguments

- object:

  An object of class `abc`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns a list with summary components.

## Details

Produces a tabular summary of an approximate Bayesian computation fit
produced by
[`ABC`](https://robustecologies.github.io/lucifer/reference/ABC.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Beaumont, M. A., Zhang, W., & Balding, D. J. (2002). Approximate
Bayesian computation in population genetics. *Genetics*, 162(4),
2025-2035.
[doi:10.1093/genetics/162.4.2025](https://doi.org/10.1093/genetics/162.4.2025)

## See also

[`plot.abc`](https://robustecologies.github.io/lucifer/reference/plot.abc.md),
[`print.abc`](https://robustecologies.github.io/lucifer/reference/print.abc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.abc
} # }
```
