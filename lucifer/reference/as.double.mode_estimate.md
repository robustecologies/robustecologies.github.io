# Coerce mode_estimate to numeric

Extracts the mode values as a numeric vector.
[`as.numeric()`](https://rdrr.io/r/base/numeric.html) dispatches to
[`as.double()`](https://rdrr.io/r/base/double.html) for primitives.

## Usage

``` r
# S3 method for class 'mode_estimate'
as.double(x, ...)
```

## Arguments

- x:

  A `mode_estimate` object.

- ...:

  Currently unused.

## Value

A numeric vector of mode values.

## Details

Implementation of `as.double.mode_estimate`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.double.mode_estimate
} # }
```
