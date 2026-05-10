# Summary method for consort objects

Extended diagnostic output including per-parameter diagnostics (MCSE,
ESS, Rhat), full condition rationale, and suggested R code.

## Usage

``` r
# S3 method for class 'consort'
summary(object, ...)
```

## Arguments

- object:

  An object of class `consort`.

- ...:

  Currently unused.

## Value

Invisibly returns `object`.

## Details

Produces a tabular summary of a Consort diagnostic report produced by
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).
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

[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`print.consort`](https://robustecologies.github.io/lucifer/reference/print.consort.md),
[`plot.consort`](https://robustecologies.github.io/lucifer/reference/plot.consort.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.consort
} # }
```
