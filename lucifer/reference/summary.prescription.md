# Summarize a prescription object

Provides the full profile details, complete scoring table, all warnings,
and interpretation guidance for a pre-fit algorithm recommendation.

## Usage

``` r
# S3 method for class 'prescription'
summary(object, n.top = 20, ...)
```

## Arguments

- object:

  An object of class `prescription`.

- n.top:

  Integer. Number of top-ranked methods to display in the scoring table.
  Default 20.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `object`.

## Details

Produces a tabular summary of a Prescribe recommendation object produced
by
[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
Yao, Y., Kennedy, L., Gabry, J., Buerkner, P.-C., & Modrak, M. (2020).
Bayesian workflow. arXiv:2011.01808.

## See also

[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
[`print.prescription`](https://robustecologies.github.io/lucifer/reference/print.prescription.md),
[`plot.prescription`](https://robustecologies.github.io/lucifer/reference/plot.prescription.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.prescription
} # }
```
