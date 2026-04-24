# Goodness-of-fit tests for cone-sample uniformity

Performs Kolmogorov-Smirnov and chi-squared goodness-of-fit tests to
check whether a sample matrix produced by
[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md)
or
[`generate_cone_samples`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
is consistent with the uniform distribution on the spherical cap of axis
\\\hat{\mu}\\ and half-angle \\\theta_0\\.

## Usage

``` r
verify_cone_uniformity(samples, mu_hat, theta0, n_bins = 20)
```

## Arguments

- samples:

  Matrix. Each column is a sample point (from `generate_cone_samples`).

- mu_hat:

  Numeric vector. Central axis direction.

- theta0:

  Numeric. Maximum planar angle of the cone.

- n_bins:

  Integer. Number of bins for histogram (default: 20).

## Value

List containing:

- ks_statistic:

  Kolmogorov-Smirnov test statistic

- ks_pvalue:

  p-value for the KS test

- chi_squared:

  Chi-squared goodness-of-fit statistic

- chi_pvalue:

  p-value for chi-squared test

- angles:

  Vector of angles between samples and \\\hat{\mu}\\

## Details

Verify uniformity of cone sampling

For each sample column, the function computes the angle \\\theta\\ to
the central axis \\\hat{\mu}\\. It then runs a Kolmogorov-Smirnov test
comparing the empirical CDF of \\\theta\\ to the theoretical
\\F\_\theta(\theta) = \Theta(\theta) / \Theta(\theta_0)\\, and a
chi-squared goodness-of-fit test on `n_bins` angular bins. Bins with
expected count below 5 are dropped before the chi-squared aggregation.
If ties in the sampled angles are detected (a finite- precision
artefact), a small jitter is applied before the KS test to avoid
spurious warnings.

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md),
[`generate_cone_samples`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md),
[`generate_hollow_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_hollow_cone_sample.md)
for the samplers under test;
[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md)
for the underlying CDF;
[`rejection_cost`](https://robustecologies.github.io/SolidAngleR/reference/rejection_cost.md)
for the expected cost of the rejection planar-angle generator.

## Examples

``` r
if (FALSE) { # \dontrun{
samples <- generate_cone_samples(10000, c(0, 0, 1), pi / 4)
fit <- verify_cone_uniformity(samples, c(0, 0, 1), pi / 4)
c(KS_p = fit$ks_pvalue, chi2_p = fit$chi_pvalue)
} # }
```
