# Skew-Laplace distribution: univariate

Density and related functions for the skew-laplace distribution:
univariate.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dslaplace(x, mu, alpha, beta, log = FALSE)

pslaplace(q, mu, alpha, beta)

qslaplace(p, mu, alpha, beta)

rslaplace(n, mu, alpha, beta)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  location parameter \\\mu\\.

- alpha:

  mixture parameter \\\alpha\\, which must be positive.

- beta:

  mixture parameter \\\beta\\, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dslaplace` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dslaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `pslaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qslaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rslaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dslaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pslaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qslaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rslaplace
} # }
```
