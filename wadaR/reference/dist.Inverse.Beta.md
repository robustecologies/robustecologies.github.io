# Inverse beta distribution

Density and random generation from the inverse beta distribution.

See Details.

See Details.

## Usage

``` r
dinvbeta(x, a, b, log = FALSE)

rinvbeta(n, a, b)
```

## Arguments

- x:

  location vector at which to evaluate density.

- a:

  scalar shape parameter \\\alpha\\.

- b:

  scalar shape parameter \\\beta\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of draws from the distribution.

## Value

`dinvbeta` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dinvbeta`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rinvbeta`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dbeta`](https://rdrr.io/r/stats/Beta.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dinvbeta
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rinvbeta
} # }
```
