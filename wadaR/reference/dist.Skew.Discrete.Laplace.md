# Skew discrete Laplace distribution: univariate

Density and related functions for the skew discrete laplace
distribution: univariate.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dsdlaplace(x, p, q, log = FALSE)

psdlaplace(x, p, q)

qsdlaplace(prob, p, q)

rsdlaplace(n, p, q)
```

## Arguments

- x:

  vector of data.

- p:

  scalar or vector of parameter \\p \in \[0,1\]\\.

- q:

  scalar or vector of parameter \\q \in \[0,1\]\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- prob:

  probability scalar or vector.

- n:

  number of observations.

## Value

`dsdlaplace` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dsdlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `psdlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qsdlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rsdlaplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dsdlaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving psdlaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qsdlaplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rsdlaplace
} # }
```
