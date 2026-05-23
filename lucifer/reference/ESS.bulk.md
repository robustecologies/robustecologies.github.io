# Bulk effective sample size

Estimates the effective sample size of rank-normalized MCMC draws,
measuring mixing in the bulk of the posterior distribution. This is the
recommended replacement for the standard
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md) when
assessing convergence, because rank-normalization makes the estimate
robust to heavy tails and infinite variance.

## Usage

``` r
ESS.bulk(x)
```

## Arguments

- x:

  a vector or matrix of posterior samples. If a matrix, each column is
  treated as a separate parameter.

## Value

A numeric vector with one bulk-ESS value per column of `x`.

## Details

The rank-normalization follows Vehtari et al. (2021): each draw is
replaced by its corresponding standard normal quantile via the Blom
(1958) fractional rank formula,

\$\$z_r = \Phi^{-1}\\\left(\frac{r - 3/8}{S - 1/4}\right)\$\$

where \\r\\ is the rank of the draw among all \\S\\ draws in that
column. The standard ESS (via AR spectral density estimation) is then
computed on the rank-normalized values.

This transformation ensures that the ESS estimate reflects mixing
behavior in the central mass of the posterior, not in the tails. For
tail mixing, see
[`ESS.tail`](https://robustecologies.github.io/lucifer/reference/ESS.tail.md).

A bulk ESS of at least 400 is generally recommended for reliable
posterior mean and median estimates. Values below 100 indicate poor
mixing in the bulk of the distribution.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and Burkner, P.-C.
(2021). "Rank-normalization, folding, and localization: An improved
R-hat for assessing convergence of MCMC (with discussion)". *Bayesian
Analysis*, 16(2), p. 667–718. doi:10.1214/20-BA1221.

## See also

[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`ESS.tail`](https://robustecologies.github.io/lucifer/reference/ESS.tail.md),
[`Rhat`](https://robustecologies.github.io/lucifer/reference/Rhat.md).

## Examples

``` r
if (FALSE) { # \dontrun{
theta <- cumsum(rnorm(2000))
ESS(theta)       # standard ESS
ESS.bulk(theta)  # bulk ESS (rank-normalized)
ESS.tail(theta)  # tail ESS (folded rank-normalized)
} # }
```
