# Summary bayes mode multi

S3 method: apply [`summary()`](https://rdrr.io/r/base/summary.html) to
objects of class `bayes_mode_multi`.

## Usage

``` r
# S3 method for class 'bayes_mode_multi'
summary(object, digits = 4, ...)
```

## Arguments

- object:

  See Details.

- digits:

  See Details.

- ...:

  See Details.

## Value

An object of class `summary.bayes_mode` with components for quantiles,
effective sample size, and diagnostic flags.

## Details

Implementation of `summary.bayes_mode_multi`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

## See also

[`plot.bayes_mode_multi`](https://robustecologies.github.io/lucifer/reference/plot.bayes_mode_multi.md),
[`print.bayes_mode_multi`](https://robustecologies.github.io/lucifer/reference/print.bayes_mode_multi.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.bayes_mode_multi
} # }
```
