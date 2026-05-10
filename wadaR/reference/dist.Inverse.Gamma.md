# Inverse gamma distribution

Density and random generation from the inverse gamma distribution.

See Details.

See Details.

## Usage

``` r
dinvgamma(x, shape = 1, scale = 1, log = FALSE)

rinvgamma(n, shape = 1, scale = 1)
```

## Arguments

- x:

  scalar location to evaluate density.

- shape:

  scalar shape parameter \\\alpha\\, default is 1.

- scale:

  scalar scale parameter \\\beta\\, default is 1.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of draws from the distribution.

## Value

`dinvgamma` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dinvgamma`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rinvgamma`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dgamma`](https://rdrr.io/r/stats/GammaDist.html),
[`dnorm`](https://rdrr.io/r/stats/Normal.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dinvgamma
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rinvgamma
} # }
```
