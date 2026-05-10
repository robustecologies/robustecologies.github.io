# Multivariate normal distribution: precision parameterization

Density and random generation for the multivariate normal distribution:
precision parameterization.

See Details.

See Details.

## Usage

``` r
dmvnp(x, mu, Omega, log = FALSE)

rmvnp(n = 1, mu = rep(0, k), Omega)
```

## Arguments

- x:

  data matrix or vector.

- mu:

  mean vector or matrix.

- Omega:

  \\K \times K\\ positive-definite precision matrix.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dmvnp` gives the density, and `rmvnp` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvnp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvnp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvnp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvnp
} # }
```
