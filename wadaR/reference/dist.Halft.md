# Half-t distribution

Density, distribution function, quantile function, and random generation
for the half-t distribution.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dhalft(x, scale = 25, nu = 1, log = FALSE)

phalft(q, scale = 25, nu = 1)

qhalft(p, scale = 25, nu = 1)

rhalft(n, scale = 25, nu = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- scale:

  scale parameter \\\alpha\\, which must be positive.

- nu:

  scalar degrees of freedom parameter \\\nu\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dhalft` gives the density or other results.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dhalft`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `phalft`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qhalft`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rhalft`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`dhalfcauchy`](https://robustecologies.github.io/lucifer/reference/dist.HalfCauchy.md),
[`dst`](https://robustecologies.github.io/lucifer/reference/dist.Student.t.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dhalft
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving phalft
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qhalft
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rhalft
} # }
```
