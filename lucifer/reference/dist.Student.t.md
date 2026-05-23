# Student t distribution: univariate

Density and related functions for the student t distribution:
univariate.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dst(x, mu = 0, sigma = 1, nu = 10, log = FALSE)

pst(q, mu = 0, sigma = 1, nu = 10, lower.tail = TRUE, log.p = FALSE)

qst(p, mu = 0, sigma = 1, nu = 10, lower.tail = TRUE, log.p = FALSE)

rst(n, mu = 0, sigma = 1, nu = 10)
```

## Arguments

- x, q:

  vector of quantiles.

- mu:

  location parameter \\\mu\\.

- sigma:

  scale parameter \\\sigma\\, which must be positive.

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

`dst` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dst`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `pst`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qst`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rst`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dst
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pst
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qst
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rst
} # }
```
