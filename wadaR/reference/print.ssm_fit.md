# Print method for SSM fit objects

Prints a concise summary of state-space model inference.

## Usage

``` r
# S3 method for class 'ssm_fit'
print(x, ...)
```

## Arguments

- x:

  An object of class `ssm_fit`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a state-space model fit
produced by
[`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Durbin, J., & Koopman, S. J. (2012). *Time Series Analysis by State
Space Methods* (2nd ed.). Oxford University Press. ISBN 9780199641178.

## See also

[`SSM`](https://robustecologies.github.io/lucifer/reference/SSM.md),
[`as.demonoid.ssm_fit`](https://robustecologies.github.io/lucifer/reference/as.demonoid.ssm_fit.md),
[`plot.ssm_fit`](https://robustecologies.github.io/lucifer/reference/plot.ssm_fit.md),
[`summary.ssm_fit`](https://robustecologies.github.io/lucifer/reference/summary.ssm_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.ssm_fit
} # }
```
