# Log-Laplace distribution: univariate symmetric

Density, distribution function, quantile function, and random generation
for the univariate symmetric log-Laplace distribution.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dllaplace(x, location = 0, scale = 1, log = FALSE)

pllaplace(q, location = 0, scale = 1)

qllaplace(p, location = 0, scale = 1)

rllaplace(n, location = 0, scale = 1)
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

`dllaplace` gives the density or other results.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dllaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `pllaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qllaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rllaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dalaplace`](https://robustecologies.github.io/lucifer/reference/dist.Asymmetric.Laplace.md),
[`dallaplace`](https://robustecologies.github.io/lucifer/reference/dist.Asymmetric.Log.Laplace.md),
[`dlaplace`](https://robustecologies.github.io/lucifer/reference/dist.Laplace.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dllaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pllaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qllaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rllaplace
} # }
```
