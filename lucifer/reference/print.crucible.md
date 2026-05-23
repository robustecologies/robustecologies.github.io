# Print method for crucible objects

Displays a compact dashboard summarizing the automated inference
pipeline results.

## Usage

``` r
# S3 method for class 'crucible'
print(x, ...)
```

## Arguments

- x:

  An object of class `crucible`.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a Crucible pipeline fit
produced by
[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
Yao, Y., Kennedy, L., Gabry, J., Buerkner, P.-C., & Modrak, M. (2020).
Bayesian workflow. arXiv:2011.01808.

## See also

[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md),
[`summary.crucible`](https://robustecologies.github.io/lucifer/reference/summary.crucible.md),
[`plot.crucible`](https://robustecologies.github.io/lucifer/reference/plot.crucible.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.crucible
} # }
```
