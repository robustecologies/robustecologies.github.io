# Summary method for crucible objects

Extended output including the Prescribe profile, per-method refinement
history, and Arena metrics.

## Usage

``` r
# S3 method for class 'crucible'
summary(object, ...)
```

## Arguments

- object:

  An object of class `crucible`.

- ...:

  Currently unused.

## Value

Invisibly returns `object`.

## Details

Produces a tabular summary of a Crucible pipeline fit produced by
[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md).
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

[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md),
[`print.crucible`](https://robustecologies.github.io/lucifer/reference/print.crucible.md),
[`plot.crucible`](https://robustecologies.github.io/lucifer/reference/plot.crucible.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.crucible
} # }
```
