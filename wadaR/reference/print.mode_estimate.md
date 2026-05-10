# Print method for mode estimates

Displays a concise summary of mode estimation results in a
bayestestR-inspired tabular format with pipe-separated columns.

## Usage

``` r
# S3 method for class 'mode_estimate'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `mode_estimate`.

- digits:

  Number of decimal digits. Default 3.

- ...:

  Currently unused.

## Value

Invisibly returns the input object.

## Details

Produces a concise one-screen console report of a mode estimate produced
by
[`Mode`](https://robustecologies.github.io/lucifer/reference/Mode.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Cross, J. L., Hoogerheide, L., Ulker, P., & van Dijk, H. K. (2024).
Sparse finite mixtures for modal inference. *Economics Letters*, 235,
111551.
[doi:10.1016/j.econlet.2024.111551](https://doi.org/10.1016/j.econlet.2024.111551)

## See also

[`Mode`](https://robustecologies.github.io/lucifer/reference/Mode.md),
[`as.double.mode_estimate`](https://robustecologies.github.io/lucifer/reference/as.double.mode_estimate.md),
[`plot.mode_estimate`](https://robustecologies.github.io/lucifer/reference/plot.mode_estimate.md),
[`summary.mode_estimate`](https://robustecologies.github.io/lucifer/reference/summary.mode_estimate.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.mode_estimate
} # }
```
