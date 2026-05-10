# Summary method for mode estimates

Provides an extended summary of mode estimation results including
modality classification, individual mode statistics, and optional
comparison across estimation methods.

## Usage

``` r
# S3 method for class 'mode_estimate'
summary(object, compare = FALSE, digits = 4, ...)
```

## Arguments

- object:

  An object of class `mode_estimate`.

- compare:

  Logical; if `TRUE`, computes modes using all available methods and
  displays a comparison. Default `FALSE`.

- digits:

  Number of decimal digits. Default 4.

- ...:

  Currently unused.

## Value

Invisibly returns a list of class `summary.mode_estimate` containing the
summary components.

## Details

Produces a tabular summary of a mode estimate produced by
[`Mode`](https://robustecologies.github.io/lucifer/reference/Mode.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Cross, J. L., Hoogerheide, L., Ulker, P., & van Dijk, H. K. (2024).
Sparse finite mixtures for modal inference. *Economics Letters*, 235,
111551.
[doi:10.1016/j.econlet.2024.111551](https://doi.org/10.1016/j.econlet.2024.111551)

## See also

[`Mode`](https://robustecologies.github.io/lucifer/reference/Mode.md),
[`as.double.mode_estimate`](https://robustecologies.github.io/lucifer/reference/as.double.mode_estimate.md),
[`plot.mode_estimate`](https://robustecologies.github.io/lucifer/reference/plot.mode_estimate.md),
[`print.mode_estimate`](https://robustecologies.github.io/lucifer/reference/print.mode_estimate.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.mode_estimate
} # }
```
