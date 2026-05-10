# Print a prescription object

Displays a concise summary of the pre-fit algorithm recommendation,
including the model fingerprint, primary recommendation with
ready-to-paste R code, and the top alternatives.

## Usage

``` r
# S3 method for class 'prescription'
print(x, ...)
```

## Arguments

- x:

  An object of class `prescription`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a Prescribe
recommendation object produced by
[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md).
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

[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
[`summary.prescription`](https://robustecologies.github.io/lucifer/reference/summary.prescription.md),
[`plot.prescription`](https://robustecologies.github.io/lucifer/reference/plot.prescription.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.prescription
} # }
```
