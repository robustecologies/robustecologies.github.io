# Multivariate Cauchy distribution: precision parameterization

Density and random generation for the multivariate cauchy distribution:
precision parameterization.

See Details.

See Details.

## Usage

``` r
dmvcp(x, mu, Omega, log = FALSE)

rmvcp(n = 1, mu, Omega)
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

`dmvcp` gives the density, and `rmvcp` generates random deviates.

See Details.

See Details.

## Details

Implementation of `dmvcp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rmvcp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dmvcp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rmvcp
} # }
```
