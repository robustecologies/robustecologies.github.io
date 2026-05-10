# Laplace distribution: precision parameterization

Density, distribution function, quantile function, and random generation
for the univariate symmetric Laplace distribution with location
parameter \\\mu\\ and precision parameter \\\tau\\.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dlaplacep(x, mu = 0, tau = 1, log = FALSE)

plaplacep(q, mu = 0, tau = 1)

qlaplacep(p, mu = 0, tau = 1)

rlaplacep(n, mu = 0, tau = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  location parameter \\\mu\\.

- tau:

  precision parameter \\\tau\\, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dlaplacep` gives the density or other results.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dlaplacep`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `plaplacep`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qlaplacep`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rlaplacep`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dlaplace`](https://robustecologies.github.io/lucifer/reference/dist.Laplace.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dlaplacep
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plaplacep
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qlaplacep
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rlaplacep
} # }
```
