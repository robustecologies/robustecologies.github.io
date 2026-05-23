# Normal-Laplace distribution: univariate asymmetric

Density and related functions for the normal-laplace distribution:
univariate asymmetric.

See Details.

See Details.

## Usage

``` r
dnormlaplace(x, mu = 0, sigma = 1, alpha = 1, beta = 1, log = FALSE)

rnormlaplace(n, mu = 0, sigma = 1, alpha = 1, beta = 1)
```

## Arguments

- x:

  vector of data.

- mu:

  location parameter \\\mu\\.

- sigma:

  scale parameter \\\sigma\\, which must be positive.

- alpha:

  shape parameter \\\alpha\\ for left-tail behavior, which must be
  positive.

- beta:

  shape parameter \\\beta\\ for right-tail behavior, which must be
  positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dnormlaplace` gives the density, and other functions provide related
computations.

See Details.

See Details.

## Details

Implementation of `dnormlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rnormlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dnormlaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rnormlaplace
} # }
```
