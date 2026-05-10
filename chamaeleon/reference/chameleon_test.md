# Statistical test for chameleon attractor behavior

Perform rigorous null hypothesis testing to determine whether a time
series exhibits chameleon attractor behavior (scale-dependent dynamical
properties). Uses surrogate-based permutation testing to establish
statistical significance.

## Usage

``` r
chameleon_test(
  x,
  n_surrogates = 199,
  surrogate_method = c("scale_shuffle", "fourier", "iaaft", "block_bootstrap"),
  test_statistics = c("D_trend", "theta_trend"),
  alpha = 0.05,
  correct_multiple = c("fisher", "bonferroni", "bh", "none"),
  robustness_checks = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  Numeric vector (time series), object of class `scale_metrics`, or
  object of class `chameleon_analysis`. If a time series is provided,
  the full analysis pipeline is run internally.

- n_surrogates:

  Integer. Number of surrogate realizations for the permutation test
  (default 199). More surrogates provide finer p-value resolution but
  increase computation time. Minimum p-value achievable is
  1/(n_surrogates + 1).

- surrogate_method:

  Character. Method for surrogate generation:

  "scale_shuffle"

  :   (Default) Shuffle frequency labels of D(f) and theta(f), breaking
      systematic scale dependence while preserving metric distributions.
      This is the primary method directly testing H0.

  "fourier"

  :   Phase randomization preserving power spectrum. Tests the linear
      null hypothesis.

  "iaaft"

  :   Iterative amplitude adjusted Fourier transform, preserving both
      power spectrum and amplitude distribution. Tests nonlinear null.

  "block_bootstrap"

  :   Block bootstrap preserving local temporal structure while
      destroying long-range dependence.

- test_statistics:

  Character vector. Statistics to compute for testing:

  "D_trend"

  :   (Recommended) Spearman correlation of D with log-frequency. Tests
      for monotonic relationship between dimension and scale.

  "theta_trend"

  :   (Recommended) Spearman correlation of theta with log-frequency.
      Tests for monotonic relationship between persistence and scale.

  "D_range"

  :   Range of mean dimension across scales. WARNING: This is
      permutation-invariant and cannot be tested with scale_shuffle
      method.

  "theta_range"

  :   Range of mean persistence across scales. WARNING:
      Permutation-invariant, not suitable for scale_shuffle.

  "D_cv"

  :   Coefficient of variation of D(f). Permutation-invariant.

  "theta_cv"

  :   Coefficient of variation of theta(f). Permutation-invariant.

  "D_max_gradient"

  :   Maximum local change in D between adjacent scales.

  "theta_max_gradient"

  :   Maximum local change in theta between adjacent scales.

  Default is `c("D_trend", "theta_trend")` which are the only statistics
  that can be meaningfully tested with the scale_shuffle method.

- alpha:

  Numeric. Significance level for hypothesis testing (default 0.05).

- correct_multiple:

  Character. Method for combining multiple test statistics:

  "fisher"

  :   (Default) Fisher's combined probability test

  "bonferroni"

  :   Bonferroni correction (conservative)

  "bh"

  :   Benjamini-Hochberg FDR correction

  "none"

  :   No correction, report individual p-values only

- robustness_checks:

  Logical. Perform additional diagnostics including embedding
  sensitivity and split-half replication (default TRUE). Increases
  computation time but improves reliability.

- parallel:

  Logical. Use parallel computation for surrogate generation (default
  TRUE). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix, sequential on Windows.

- n_cores:

  Integer. Number of CPU cores for parallel computation (default is
  `parallel::detectCores() - 1`).

- verbose:

  Logical. Print progress information (default FALSE).

- ...:

  Additional arguments passed to internal functions (e.g., embedding
  parameters m, tau if x is a raw time series).

## Value

Object of class `"chameleon_test"` containing:

- observed:

  Named list of observed test statistics.

- null_distributions:

  Matrix with surrogate test statistics (n_surrogates x n_statistics).

- p_values:

  Named numeric vector of p-values per statistic.

- p_combined:

  Numeric. Combined p-value using specified correction method.

- effect_sizes:

  Data frame with Cohen's d, percentile rank, and SNR for each
  statistic.

- conclusion:

  Character. Either "chameleon" or "non-chameleon".

- confidence:

  Character. Confidence level: "high", "moderate", or "low".

- diagnostics:

  List with robustness check results (if performed).

- scale_metrics:

  The scale_metrics object used for testing.

- params:

  List of test parameters.

S3 methods available: [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Details

**Null hypothesis:**

H0: D(f) and theta(f) are constant across frequency scales (no chameleon
behavior).

H1: D(f) and/or theta(f) vary systematically across frequency scales.

**Testing procedure:**

1.  Compute observed test statistics from the scale-dependent metrics.

2.  Generate n_surrogates realizations under the null hypothesis.

3.  Compute the same test statistics for each surrogate.

4.  Calculate p-values as the proportion of surrogate statistics \>=
    observed.

5.  Apply multiple testing correction if multiple statistics are used.

6.  Compute effect sizes to quantify practical significance.

**P-value computation:**

Uses the conservative formula from Phipson & Smyth (2010): \$\$p =
\frac{\sum\_{i=1}^{B} I(\|T_i^\*\| \geq \|T\_{obs}\|) + 1}{B + 1}\$\$

This avoids p = 0 and provides valid permutation p-values.

**Effect size interpretation (Cohen's d):**

- \|d\| \< 0.5: Small effect

- 0.5 \<= \|d\| \< 0.8: Medium effect

- \|d\| \>= 0.8: Large effect

- \|d\| \>= 2.0: Very large effect

**Robustness checks:**

- **Embedding sensitivity**: Re-test with m +/- 1, tau +/- 2

- **Split-half**: Test first and second halves separately

- Warnings issued if results are sensitive to parameter choices

## References

Phipson B, Smyth GK (2010). Permutation p-values should never be zero:
calculating exact p-values when permutations are randomly drawn.
Statistical Applications in Genetics and Molecular Biology 9:Article 39.
[doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585)

Theiler J, Eubank S, Longtin A, Galdrikian B, Farmer JD (1992). Testing
for nonlinearity in time series: the method of surrogate data. Physica D
58:77-94.
[doi:10.1016/0167-2789(92)90102-S](https://doi.org/10.1016/0167-2789%2892%2990102-S)

Schreiber T, Schmitz A (2000). Surrogate time series. Physica D
142:346-382.
[doi:10.1016/S0167-2789(00)00043-9](https://doi.org/10.1016/S0167-2789%2800%2900043-9)

Fisher RA (1932). Statistical methods for research workers. Oliver and
Boyd, Edinburgh. ISBN: 978-0-05-002170-5.

## See also

[`generate_surrogates`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md)
for surrogate generation methods,
[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
for the complete analysis workflow,
[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
for computing scale metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
# Test a time series for chameleon behavior
set.seed(42)
t <- seq(0, 100, by = 0.01)
x <- sin(2*pi*0.1*t) + 0.5*sin(2*pi*t) + 0.2*rnorm(length(t))

# Run statistical test
test_result <- chameleon_test(x, n_surrogates = 99, verbose = TRUE)
print(test_result)
plot(test_result)

# Test existing scale_metrics object
analysis <- chameleon_analysis(x)
test_result <- chameleon_test(analysis$scale_metrics)
} # }
```
