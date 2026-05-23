# (Scaled) inverse chi-squared distribution

Density and random generation for the (scaled) inverse chi-squared
distribution.

See Details.

See Details.

## Usage

``` r
dinvchisq(x, df, scale = 1/df, log = FALSE)

rinvchisq(n, df, scale = 1/df)
```

## Arguments

- x:

  vector of quantiles.

- df:

  degrees of freedom parameter \\\nu\\.

- scale:

  scale parameter \\\lambda\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations.

## Value

`dinvchisq` gives the density or other results.

See Details.

See Details.

## Details

Implementation of `dinvchisq`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rinvchisq`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dchisq`](https://rdrr.io/r/stats/Chisquare.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dinvchisq
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rinvchisq
} # }
```
