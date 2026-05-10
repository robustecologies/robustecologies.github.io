# Print method for node_fit objects

S3 method: apply [`print()`](https://rdrr.io/r/base/print.html) to
objects of class `node_fit`.

## Usage

``` r
# S3 method for class 'node_fit'
print(x, ...)
```

## Arguments

- x:

  An object of class `node_fit`.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a neural ODE fit
produced by
[`NODE`](https://robustecologies.github.io/lucifer/reference/NODE.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Bonnaffe, W., Sheldon, B. C., & Coulson, T. (2023). Neural ordinary
differential equations for ecological and evolutionary time series
analysis. *Methods in Ecology and Evolution*, 12(7), 1301-1315.
[doi:10.1111/2041-210X.13606](https://doi.org/10.1111/2041-210X.13606)

## See also

[`NODE`](https://robustecologies.github.io/lucifer/reference/NODE.md),
[`plot.node_fit`](https://robustecologies.github.io/lucifer/reference/plot.node_fit.md),
[`predict.node_fit`](https://robustecologies.github.io/lucifer/reference/NODE_predict.md),
[`summary.node_fit`](https://robustecologies.github.io/lucifer/reference/summary.node_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.node_fit
} # }
```
