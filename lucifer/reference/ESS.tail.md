# Tail effective sample size

Estimates the effective sample size for the tails of the posterior
distribution by folding draws around their median, rank-normalizing, and
computing ESS. This detects poor mixing in the tails that both the
standard
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md) and
[`ESS.bulk`](https://robustecologies.github.io/lucifer/reference/ESS.bulk.md)
may miss.

## Usage

``` r
ESS.tail(x)
```

## Arguments

- x:

  a vector or matrix of posterior samples. If a matrix, each column is
  treated as a separate parameter.

## Value

A numeric vector with one tail-ESS value per column of `x`.

## Details

The folding transformation of Vehtari et al. (2021) maps each draw to
its absolute deviation from the median:

\$\$x^{\mathrm{fold}}\_i = \|x_i - \mathrm{median}(x)\|\$\$

This amplifies differences in the tails, making ESS sensitive to tail
mixing. The folded draws are then rank-normalized using the same Blom
formula as
[`ESS.bulk`](https://robustecologies.github.io/lucifer/reference/ESS.bulk.md)
and ESS is computed on the result.

When the chain mixes poorly in the tails (e.g., slow exploration of
extreme values), the folded draws will have high autocorrelation and the
tail ESS will be small. This is critical for reliable estimation of
posterior quantiles and credible intervals.

A tail ESS of at least 400 is generally recommended for reliable
posterior interval estimates. Values below 100 suggest that tail
quantiles (e.g., 2.5% and 97.5%) are unreliable.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and Burkner, P.-C.
(2021). "Rank-normalization, folding, and localization: An improved
R-hat for assessing convergence of MCMC (with discussion)". *Bayesian
Analysis*, 16(2), p. 667–718. doi:10.1214/20-BA1221.

## See also

[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`ESS.bulk`](https://robustecologies.github.io/lucifer/reference/ESS.bulk.md),
[`Rhat`](https://robustecologies.github.io/lucifer/reference/Rhat.md).

## Examples

``` r
if (FALSE) { # \dontrun{
theta <- cumsum(rnorm(2000))
ESS(theta)       # standard ESS
ESS.bulk(theta)  # bulk ESS
ESS.tail(theta)  # tail ESS (folded rank-normalized)
} # }
```
