# Dirichlet distribution

Density and random generation from the Dirichlet distribution.

See Details.

See Details.

## Usage

``` r
ddirichlet(x, alpha, log = FALSE)

rdirichlet(n, alpha)
```

## Arguments

- x:

  vector or matrix containing one random deviate per row, each summing
  to 1.

- alpha:

  vector or matrix of shape parameters.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of random deviates to generate.

## Value

`ddirichlet` gives the density, and `rdirichlet` generates random
deviates.

See Details.

See Details.

## Details

The Dirichlet distribution is the multivariate generalization of the
univariate beta distribution. It is the conjugate prior distribution for
the parameters of the categorical and multinomial distributions.

Implementation of `ddirichlet`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rdirichlet`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dbeta`](https://rdrr.io/r/stats/Beta.html),
[`dcat`](https://robustecologies.github.io/lucifer/reference/dist.Categorical.md),
[`dmvpolya`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Polya.md),
[`dmultinom`](https://rdrr.io/r/stats/Multinom.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- ddirichlet(c(.1,.3,.6), c(1,1,1))
x <- rdirichlet(10, c(1,1,1))
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving ddirichlet
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rdirichlet
} # }
```
