# Power exponential distribution: univariate symmetric

Density and related functions for the power exponential distribution:
univariate symmetric.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dpe(x, mu = 0, sigma = 1, kappa = 2, log = FALSE)

ppe(q, mu = 0, sigma = 1, kappa = 2, lower.tail = TRUE, log.p = FALSE)

qpe(p, mu = 0, sigma = 1, kappa = 2, lower.tail = TRUE, log.p = FALSE)

rpe(n, mu = 0, sigma = 1, kappa = 2)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  location parameter \\\mu\\.

- sigma:

  scale parameter \\\sigma\\, which must be positive.

- kappa:

  kurtosis parameter \\\kappa\\, which must be positive.

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

`dpe` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dpe`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `ppe`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qpe`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rpe`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dpe
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving ppe
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qpe
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rpe
} # }
```
