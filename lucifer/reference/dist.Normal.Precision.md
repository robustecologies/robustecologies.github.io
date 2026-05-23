# Normal distribution: precision parameterization

Density and related functions for the normal distribution: precision
parameterization.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dnormp(x, mean = 0, prec = 1, log = FALSE)

pnormp(q, mean = 0, prec = 1, lower.tail = TRUE, log.p = FALSE)

qnormp(p, mean = 0, prec = 1, lower.tail = TRUE, log.p = FALSE)

rnormp(n, mean = 0, prec = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- mean:

  mean parameter.

- prec:

  precision parameter, which must be positive.

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

`dnormp` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dnormp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `pnormp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qnormp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rnormp`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dnormp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pnormp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qnormp
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rnormp
} # }
```
