# Generalized Pareto distribution

Density and random generation for the generalized Pareto distribution.

See Details.

See Details.

## Usage

``` r
dgpd(x, mu, sigma, xi, log = FALSE)

rgpd(n, mu, sigma, xi)
```

## Arguments

- x:

  vector of data.

- mu:

  scalar or vector location parameter \\\mu\\.

- sigma:

  positive-only scalar or vector scale parameter \\\sigma\\.

- xi:

  scalar or vector shape parameter \\\xi\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  positive scalar integer, the number of observations to generate.

## Value

`dgpd` gives the density, and `rgpd` generates random deviates.

See Details.

See Details.

## Details

The generalized Pareto distribution (GPD) is a more flexible extension
of the Pareto distribution. It is equivalent to the exponential
distribution when \\\mu = 0\\ and \\\xi = 0\\.

Implementation of `dgpd`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rgpd`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Pickands J. (1975). "Statistical Inference Using Extreme Order
Statistics". *The Annals of Statistics*, 3, p. 119–131.

## See also

[`dpareto`](https://robustecologies.github.io/lucifer/reference/dist.Pareto.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- dgpd(0, 0, 1, 0, log=TRUE)
x <- rgpd(10, 0, 1, 0)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dgpd
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rgpd
} # }
```
