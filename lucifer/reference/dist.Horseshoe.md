# Horseshoe distribution

Density and random generation from the horseshoe distribution.

See Details.

See Details.

## Usage

``` r
dhs(x, lambda, tau, log = FALSE)

rhs(n, lambda, tau)
```

## Arguments

- x:

  location vector at which to evaluate density.

- lambda:

  positive-only local parameter vector \\\lambda\\.

- tau:

  positive-only global parameter scalar \\\tau\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of draws from the distribution.

## Value

`dhs` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dhs`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rhs`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`dlaplace`](https://robustecologies.github.io/lucifer/reference/dist.Laplace.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dhs
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rhs
} # }
```
