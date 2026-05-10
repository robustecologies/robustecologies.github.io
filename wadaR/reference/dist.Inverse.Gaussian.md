# Inverse Gaussian distribution

Density and random generation from the inverse Gaussian distribution.

See Details.

See Details.

## Usage

``` r
dinvgaussian(x, mu, lambda, log = FALSE)

rinvgaussian(n, mu, lambda)
```

## Arguments

- x:

  scalar location to evaluate density.

- mu:

  mean parameter \\\mu\\, which must be positive.

- lambda:

  inverse-variance parameter \\\lambda\\, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of draws from the distribution.

## Value

`dinvgaussian` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dinvgaussian`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rinvgaussian`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dnorm`](https://rdrr.io/r/stats/Normal.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dinvgaussian
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rinvgaussian
} # }
```
