# Pareto distribution

Density and related functions for the pareto distribution.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dpareto(x, alpha, log = FALSE)

ppareto(q, alpha)

qpareto(p, alpha)

rpareto(n, alpha)
```

## Arguments

- x, q:

  vector of quantiles.

- alpha:

  shape parameter \\\alpha\\, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dpareto` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dpareto`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `ppareto`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qpareto`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rpareto`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dpareto
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving ppareto
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qpareto
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rpareto
} # }
```
