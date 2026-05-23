# Simplex distribution

Density and random generation for the simplex distribution on \\(0,
1)\\.

See Details.

See Details.

## Usage

``` r
dsimplex(x, mu = 0.5, sigma = 1, log = FALSE)

rsimplex(n, mu = 0.5, sigma = 1)
```

## Arguments

- x:

  vector of quantiles in \\(0, 1)\\.

- mu:

  vector of mean parameters in \\(0, 1)\\.

- sigma:

  vector of dispersion parameters, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations, which must be a positive integer of length 1.

## Value

`dsimplex` gives the density and `rsimplex` generates random deviates.

See Details.

See Details.

## Details

The simplex distribution (Barndorff-Nielsen and Jorgensen, 1991) is a
two-parameter distribution on \\(0, 1)\\ that provides an alternative to
the beta distribution for proportional data. Its density is \$\$f(y;
\mu, \sigma) = \[2\pi\sigma^2 y^3 (1-y)^3\]^{-1/2} \exp\left(-\frac{d(y;
\mu)}{2\sigma^2}\right)\$\$ where \\d(y; \mu) = (y - \mu)^2 /
\[y(1-y)\mu^2(1-\mu)^2\]\\.

Implementation of `dsimplex`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rsimplex`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Barndorff-Nielsen, O.E. and Jorgensen, B. (1991). "Some Parametric
Models on the Simplex". *Journal of Multivariate Analysis*, 39(1), p.
106–116.

## See also

[`dbeta`](https://rdrr.io/r/stats/Beta.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- seq(0.01, 0.99, length.out = 200)
d <- dsimplex(x, mu = 0.3, sigma = 1)
plot(x, d, type = "l")
r <- rsimplex(1000, mu = 0.5, sigma = 0.5)
hist(r, breaks = 50, freq = FALSE)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dsimplex
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rsimplex
} # }
```
