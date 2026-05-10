# Log-normal distribution: precision parameterization

Density, distribution function, quantile function, and random generation
for the univariate log-normal distribution with mean \\\mu\\ and
precision \\\tau\\.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dlnormp(x, mu, tau = NULL, var = NULL, log = FALSE)

plnormp(q, mu, tau, lower.tail = TRUE, log.p = FALSE)

qlnormp(p, mu, tau, lower.tail = TRUE, log.p = FALSE)

rlnormp(n, mu, tau = NULL, var = NULL)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  mean parameter \\\mu\\.

- tau:

  precision parameter \\\tau\\, which must be positive. Cannot be used
  together with `var`.

- var:

  variance parameter, which must be positive. Cannot be used together
  with `tau`.

- log, log.p:

  logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\Pr\[X \le x\]\\.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dlnormp` gives the density or other results.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dlnormp`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `plnormp`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qlnormp`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rlnormp`. Refer to the package vignettes and the
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
## see package vignettes for full examples involving dlnormp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plnormp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qlnormp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rlnormp
} # }
```
