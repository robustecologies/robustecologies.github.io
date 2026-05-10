# Normal distribution: variance parameterization

Density and related functions for the normal distribution: variance
parameterization.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dnormv(x, mean = 0, var = 1, log = FALSE)

pnormv(q, mean = 0, var = 1, lower.tail = TRUE, log.p = FALSE)

qnormv(p, mean = 0, var = 1, lower.tail = TRUE, log.p = FALSE)

rnormv(n, mean = 0, var = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- mean:

  mean parameter.

- var:

  variance parameter, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\Pr\[X \le x\]\\.

- log.p:

  logical; if `TRUE`, probabilities are given as log(p).

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dnormv` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dnormv`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `pnormv`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qnormv`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rnormv`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dnormv
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pnormv
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qnormv
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rnormv
} # }
```
