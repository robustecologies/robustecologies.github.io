# Normal distribution mixture

Density and related functions for the normal distribution mixture.

See Details.

See Details.

See Details.

## Usage

``` r
dnormm(x, p, mu, sigma, log = FALSE)

pnormm(q, p, mu, sigma, lower.tail = TRUE, log.p = FALSE)

rnormm(n, p, mu, sigma)
```

## Arguments

- x, q:

  vector of quantiles.

- p:

  vector of length \\M\\ of probabilities for components.

- mu:

  vector of length \\M\\ of means.

- sigma:

  vector of length \\M\\ of positive standard deviations.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\Pr\[X \le x\]\\.

- log.p:

  logical; if `TRUE`, probabilities are given as log(p).

- n:

  number of observations.

## Value

`dnormm` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

## Details

Implementation of `dnormm`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `pnormm`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rnormm`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dnormm
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pnormm
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rnormm
} # }
```
