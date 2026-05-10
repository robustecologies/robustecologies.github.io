# Half-normal distribution

Density, distribution function, quantile function, and random generation
for the half-normal distribution.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dhalfnorm(x, scale = sqrt(pi/2), log = FALSE)

phalfnorm(q, scale = sqrt(pi/2), lower.tail = TRUE, log.p = FALSE)

qhalfnorm(p, scale = sqrt(pi/2), lower.tail = TRUE, log.p = FALSE)

rhalfnorm(n, scale = sqrt(pi/2))
```

## Arguments

- x, q:

  vector of quantiles.

- scale:

  scale parameter \\\sigma\\, which must be positive.

- log, log.p:

  logical; if `TRUE`, the logarithm of the density or result is
  returned.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\Pr\[X \le x\]\\.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dhalfnorm` gives the density or other results.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dhalfnorm`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `phalfnorm`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qhalfnorm`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rhalfnorm`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dnorm`](https://rdrr.io/r/stats/Normal.html),
[`dnormp`](https://robustecologies.github.io/lucifer/reference/dist.Normal.Precision.md),
[`dnormv`](https://robustecologies.github.io/lucifer/reference/dist.Normal.Variance.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dhalfnorm
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving phalfnorm
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qhalfnorm
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rhalfnorm
} # }
```
