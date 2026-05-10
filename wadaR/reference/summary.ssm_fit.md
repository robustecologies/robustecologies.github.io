# Summary method for SSM fit objects

Provides extended diagnostics for state-space model inference.

## Usage

``` r
# S3 method for class 'ssm_fit'
summary(object, ...)
```

## Arguments

- object:

  An object of class `ssm_fit`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns a list of diagnostics.

## Details

Produces a tabular summary of a state-space model fit produced by
[`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Durbin, J., & Koopman, S. J. (2012). *Time Series Analysis by State
Space Methods* (2nd ed.). Oxford University Press. ISBN 9780199641178.

## See also

[`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md),
[`as.demonoid.ssm_fit`](https://robustecologies.github.io/lucifer/reference/as.demonoid.ssm_fit.md),
[`plot.ssm_fit`](https://robustecologies.github.io/lucifer/reference/plot.ssm_fit.md),
[`print.ssm_fit`](https://robustecologies.github.io/lucifer/reference/print.ssm_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.ssm_fit
} # }
```
