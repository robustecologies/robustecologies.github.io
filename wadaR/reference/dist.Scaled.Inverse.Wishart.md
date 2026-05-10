# Scaled inverse Wishart distribution

Density and related functions for the scaled inverse wishart
distribution.

See Details.

See Details.

## Usage

``` r
dsiw(Q, nu, S, zeta, mu, delta, log = FALSE)

rsiw(nu, S, mu, delta)
```

## Arguments

- Q:

  symmetric, positive-definite \\k \times k\\ matrix.

- nu:

  scalar degrees of freedom \\\nu\\.

- S:

  symmetric, positive-semidefinite \\k \times k\\ scale matrix.

- zeta:

  positive-only vector of length \\k\\ of auxiliary scale parameters.

- mu:

  vector of length \\k\\ of location hyperparameters.

- delta:

  positive-only vector of length \\k\\ of scale hyperparameters.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dsiw` gives the density, and other functions provide related
computations.

See Details.

See Details.

## Details

Implementation of `dsiw`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rsiw`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dsiw
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rsiw
} # }
```
