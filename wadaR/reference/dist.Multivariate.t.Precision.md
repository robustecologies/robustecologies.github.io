# Multivariate t distribution: precision parameterization

Density and random generation for the multivariate t distribution:
precision parameterization.

See Details.

See Details.

## Usage

``` r
dmvtp(x, mu, Omega, nu = Inf, log = FALSE)

rmvtp(n = 1, mu, Omega, nu = Inf)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- Omega:

  \\K \times K\\ positive-definite precision matrix.

- nu:

  positive degrees of freedom parameter.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvtp` gives the density, and `rmvtp` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvtp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvtp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvtp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvtp
} # }
```
