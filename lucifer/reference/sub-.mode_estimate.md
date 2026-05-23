# Subscript mode_estimate

Extracts mode values by index for backward compatibility. Ensures
`Mode(x)[1]` still returns the scalar mode value.

## Usage

``` r
# S3 method for class 'mode_estimate'
x[i, ...]
```

## Arguments

- x:

  A `mode_estimate` object.

- i:

  Index.

- ...:

  Currently unused.

## Value

The mode value(s) at the specified index.

## Details

Implementation of `[.mode_estimate`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving [.mode_estimate
} # }
```
