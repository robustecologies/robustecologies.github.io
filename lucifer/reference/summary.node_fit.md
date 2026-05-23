# Summary method for node_fit objects

S3 method: apply [`summary()`](https://rdrr.io/r/base/summary.html) to
objects of class `node_fit`.

## Usage

``` r
# S3 method for class 'node_fit'
summary(object, ...)
```

## Arguments

- object:

  An object of class `node_fit`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `summary.node_fit` (invisibly), printed with extended
diagnostics.

## Details

Produces a tabular summary of a neural ODE fit produced by
[`NODE`](https://robustecologies.github.io/lucifer/reference/NODE.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Bonnaffe, W., Sheldon, B. C., & Coulson, T. (2023). Neural ordinary
differential equations for ecological and evolutionary time series
analysis. *Methods in Ecology and Evolution*, 12(7), 1301-1315.
[doi:10.1111/2041-210X.13606](https://doi.org/10.1111/2041-210X.13606)

## See also

[`NODE`](https://robustecologies.github.io/lucifer/reference/NODE.md),
[`plot.node_fit`](https://robustecologies.github.io/lucifer/reference/plot.node_fit.md),
[`predict.node_fit`](https://robustecologies.github.io/lucifer/reference/NODE_predict.md),
[`print.node_fit`](https://robustecologies.github.io/lucifer/reference/print.node_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.node_fit
} # }
```
