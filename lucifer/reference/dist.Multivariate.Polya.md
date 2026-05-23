# Multivariate Polya distribution

Density and random generation for the multivariate polya distribution.

See Details.

See Details.

## Usage

``` r
dmvpolya(x, alpha, log = FALSE)

rmvpolya(n = 1, alpha)
```

## Arguments

- x:

  vector of counts.

- alpha:

  vector of concentration parameters.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvpolya` gives the density, and `rmvpolya` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvpolya`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rmvpolya`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvpolya
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvpolya
} # }
```
