# Print method for multi-parameter Bayesian mode inference

S3 method: apply [`print()`](https://rdrr.io/r/base/print.html) to
objects of class `bayes_mode_multi`.

## Usage

``` r
# S3 method for class 'bayes_mode_multi'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `bayes_mode_multi`.

- digits:

  Number of decimal digits. Default 3.

- ...:

  Currently unused.

## Value

Invisibly returns `x`. The side effect is the printed report.

## Details

Implementation of `print.bayes_mode_multi`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

## See also

[`plot.bayes_mode_multi`](https://robustecologies.github.io/lucifer/reference/plot.bayes_mode_multi.md),
[`summary.bayes_mode_multi`](https://robustecologies.github.io/lucifer/reference/summary.bayes_mode_multi.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.bayes_mode_multi
} # }
```
