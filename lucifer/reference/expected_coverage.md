# Expected coverage diagnostic for SBI

Checks whether the credible intervals from the posterior have the
correct coverage. For each of `n_sims` synthetic datasets, the method
draws posterior samples and checks whether the true parameter falls
within HPD intervals at multiple credibility levels.

## Usage

``` r
expected_coverage(
  sbi_fit,
  simulator,
  prior,
  n_sims = 300L,
  n_posterior_samples = 500L,
  levels = seq(0.1, 0.9, by = 0.1),
  summary_fn = NULL,
  verbose = TRUE
)
```

## Arguments

- sbi_fit:

  An object of class `sbi`.

- simulator:

  The simulator function.

- prior:

  The prior function.

- n_sims:

  Integer: number of calibration datasets. Default 300.

- n_posterior_samples:

  Integer: posterior samples per dataset. Default 500.

- levels:

  Numeric vector of credibility levels to check. Default
  `seq(0.1, 0.9, by = 0.1)`.

- summary_fn:

  Optional summary function.

- verbose:

  Logical: print progress. Default `TRUE`.

## Value

A list of class `sbi_coverage` containing:

- coverage:

  Matrix (length(levels) x d) of empirical coverage fractions.

- levels:

  The credibility levels tested.

- n_sims:

  Number of calibration datasets.

## Details

Implementation of `expected_coverage`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## References

Deistler, M. et al. (2025). Simulation-based inference: a practical
guide. *arXiv:2411.17337*.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- SBI(simulator, prior, x_obs, n_simulations = 5000)
cov <- expected_coverage(fit, simulator, prior, n_sims = 200)
plot(fit, type = "coverage", coverage = cov)
} # }
```
