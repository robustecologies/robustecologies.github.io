# Continuous relaxation of a Markov random field distribution

Density and random generation from the continuous relaxation of a Markov
random field (MRF) distribution.

See Details.

See Details.

## Usage

``` r
dcrmrf(x, alpha, Omega, log = FALSE)

rcrmrf(n = 1, alpha, Omega)
```

## Arguments

- x:

  vector of length \\k\\.

- alpha:

  vector of length \\k\\ of shape parameters.

- Omega:

  \\k \times k\\ precision matrix \\\Omega\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of random deviates to generate.

## Value

`dcrmrf` gives the density, and `rcrmrf` generates random deviates.

See Details.

See Details.

## Details

Zhang et al. (2012) introduced a generalized form of the Gaussian
integral trick from statistical physics to transform a discrete variable
so that it may be estimated with continuous variables. An auxiliary
Gaussian variable is added to a discrete Markov random field (MRF) so
that discrete dependencies cancel out.

Implementation of `dcrmrf`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rcrmrf`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Zhang, Y., Ghahramani, Z., Storkey, A.J., and Sutton, C.A. (2012).
"Continuous Relaxations for Discrete Hamiltonian Monte Carlo". *Advances
in Neural Information Processing Systems*, 25, p. 3203–3211.

## See also

[`dmvn`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Normal.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- dcrmrf(rnorm(5), rnorm(5), diag(5))
x <- rcrmrf(10, rnorm(5), diag(5))
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dcrmrf
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rcrmrf
} # }
```
