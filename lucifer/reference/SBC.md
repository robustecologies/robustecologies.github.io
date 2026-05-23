# Simulation-based calibration for SBI

Validates the calibration of a trained SBI model by checking whether
posterior samples are consistent with the prior predictive distribution.
For each of `n_sims` synthetic datasets generated from the prior, the
method draws posterior samples using the trained network and computes
the rank of the true parameter within those samples. If the posterior is
well-calibrated, these ranks should be uniformly distributed.

## Usage

``` r
SBC(
  sbi_fit,
  simulator,
  prior,
  n_sims = 300L,
  n_posterior_samples = 200L,
  summary_fn = NULL,
  verbose = TRUE
)
```

## Arguments

- sbi_fit:

  An object of class `sbi`, as returned by
  [`SBI`](https://robustecologies.github.io/lucifer/reference/SBI.md).

- simulator:

  The same simulator function used in
  [`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md).

- prior:

  The same prior function used in
  [`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md).

- n_sims:

  Integer: number of synthetic calibration datasets to generate. Default
  300.

- n_posterior_samples:

  Integer: number of posterior samples to draw per calibration dataset.
  Default 200.

- summary_fn:

  Optional summary function (same as used in
  [`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)).

- verbose:

  Logical: print progress. Default `TRUE`.

## Value

A list of class `sbi_sbc` containing:

- ranks:

  Matrix (n_sims x d) of rank statistics.

- ks_pvalues:

  Kolmogorov-Smirnov test p-values against Uniform(0,
  n_posterior_samples) for each parameter.

- n_sims:

  Number of calibration datasets.

- n_posterior_samples:

  Number of posterior samples per dataset.

## Details

For each calibration iteration \\i\\:

1.  Draw \\\theta^\* \sim p(\theta)\\ from the prior.

2.  Simulate \\x^\* \sim p(x \mid \theta^\*)\\.

3.  Draw \\L\\ posterior samples \\\theta_1, \ldots, \theta_L \sim
    q(\theta \mid x^\*)\\.

4.  Compute the rank \\r_i = \sum\_{l=1}^{L} \mathbf{1}\[\theta_l \<
    \theta^\*\]\\.

Under correct calibration, \\r_i \sim \text{Uniform}(0, L)\\ for each
parameter dimension. The returned ranks can be visualized as histograms
via `plot(sbi_fit, type = "sbc")`.

## References

Talts, S., Betancourt, M., Simpson, D., Vehtari, A. & Gelman, A. (2018).
Validating Bayesian inference algorithms with simulation-based
calibration. *arXiv:1804.06788*.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
simulator <- function(theta) rnorm(10, mean = theta[1], sd = 1)
prior <- function() rnorm(1, mean = 0, sd = 5)
x_obs <- rnorm(10, mean = 2, sd = 1)

fit <- SBI(simulator, prior, x_obs, n_simulations = 5000)
sbc <- SBC(fit, simulator, prior, n_sims = 200)
plot(fit, type = "sbc", sbc = sbc)
} # }
```
