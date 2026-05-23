# Mixture of Laplace distributions

Density, cumulative distribution, and random generation for the mixture
of univariate Laplace distributions.

See Details.

See Details.

See Details.

## Usage

``` r
dlaplacem(x, p, location, scale, log = FALSE)

plaplacem(q, p, location, scale)

rlaplacem(n, p, location, scale)
```

## Arguments

- x, q:

  vector of values at which the density will be evaluated.

- p:

  vector of length \\M\\ of probabilities for \\M\\ components, summing
  to one.

- location:

  vector of length \\M\\ of location parameters.

- scale:

  vector of length \\M\\ of positive scale parameters.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dlaplacem` gives the density or other results.

See Details.

See Details.

See Details.

## Details

Implementation of `dlaplacem`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `plaplacem`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rlaplacem`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`ddirichlet`](https://robustecologies.github.io/lucifer/reference/dist.Dirichlet.md),
[`dlaplace`](https://robustecologies.github.io/lucifer/reference/dist.Laplace.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dlaplacem
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plaplacem
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rlaplacem
} # }
```
