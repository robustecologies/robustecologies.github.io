# Yang-Berger distribution

Density and related functions for the yang-berger distribution.

See Details.

See Details.

## Usage

``` r
dyangberger(x, log = FALSE)

dyangbergerc(x, log = FALSE)
```

## Arguments

- x:

  \\k \times k\\ positive-definite covariance or precision matrix (or
  Cholesky factor for `dyangbergerc`).

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dyangberger` gives the density, and other functions provide related
computations.

See Details.

See Details.

## Details

Implementation of `dyangberger`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `dyangbergerc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dyangberger
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dyangbergerc
} # }
```
