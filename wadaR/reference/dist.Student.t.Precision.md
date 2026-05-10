# Student t distribution: precision parameterization

Density and related functions for the student t distribution: precision
parameterization.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dstp(x, mu = 0, tau = 1, nu = 10, log = FALSE)

pstp(q, mu = 0, tau = 1, nu = 10, lower.tail = TRUE, log.p = FALSE)

qstp(p, mu = 0, tau = 1, nu = 10, lower.tail = TRUE, log.p = FALSE)

rstp(n, mu = 0, tau = 1, nu = 10)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  location parameter \\\mu\\.

- tau:

  precision parameter \\\tau\\, which must be positive.

- nu:

  degrees of freedom parameter \\\nu\\, which must be positive.

- log, log.p:

  logical; if `TRUE`, the logarithm of the density or probability is
  returned.

- lower.tail:

  logical; if `TRUE`, probabilities are \\Pr\[X \le x\]\\.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dstp` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dstp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `pstp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qstp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rstp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dstp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pstp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qstp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rstp
} # }
```
