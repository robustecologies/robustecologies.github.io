# Laplace distribution: univariate symmetric

Density, distribution function, quantile function, and random generation
for the univariate symmetric Laplace distribution with location
parameter \\\mu\\ and scale parameter \\\lambda\\.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dlaplace(x, location = 0, scale = 1, log = FALSE)

plaplace(q, location = 0, scale = 1)

qlaplace(p, location = 0, scale = 1)

rlaplace(n, location = 0, scale = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- location:

  location parameter \\\mu\\.

- scale:

  scale parameter \\\lambda\\, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dlaplace` gives the density or other results.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `plaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dalaplace`](https://robustecologies.github.io/lucifer/reference/dist.Asymmetric.Laplace.md),
[`dallaplace`](https://robustecologies.github.io/lucifer/reference/dist.Asymmetric.Log.Laplace.md),
[`dlaplacep`](https://robustecologies.github.io/lucifer/reference/dist.Laplace.Precision.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dlaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qlaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rlaplace
} # }
```
