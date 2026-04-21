# Verify uniformity of cone sampling

Performs statistical tests to verify the uniformity of samples generated
on a spherical cap.

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

This function computes the angle \\\theta\\ between each sample and the
cone axis, then performs:

1.  Kolmogorov-Smirnov test comparing the empirical CDF to the
    theoretical CDF

2.  Chi-squared goodness-of-fit test on binned angle frequencies

The theoretical CDF for angles is: \$\$F\_\theta(\theta) =
\frac{\Theta(\theta)}{\Theta(\theta_0)}\$\$

If ties are detected in the sampled angles (finite-precision effects), a
small jitter is applied before the KS test to avoid spurious warnings.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate samples and verify uniformity
samples <- generate_cone_samples(10000, c(0, 0, 1), pi/4)
results <- verify_cone_uniformity(samples, c(0, 0, 1), pi/4)

cat("KS statistic:", results$ks_statistic, "\n")
cat("KS p-value:", results$ks_pvalue, "\n")

# Plot histogram
hist(results$angles, breaks = 30, probability = TRUE,
     main = "Distribution of angles from cone axis",
     xlab = expression(theta ~ "(radians)"))
} # }
```
