# Tests of accuracy with random points (TARP)

A more powerful calibration diagnostic than SBC, based on checking
whether the true parameter falls within posterior credible regions
defined by distances to random reference points. Detects both global and
local miscalibration.

## Usage

``` r
TARP(
  sbi_fit,
  simulator,
  prior,
  n_sims = 300L,
  n_posterior_samples = 500L,
  n_references = 100L,
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

  Integer: number of synthetic datasets. Default 300.

- n_posterior_samples:

  Integer: posterior samples per dataset. Default 500.

- n_references:

  Integer: number of random reference points. Default 100.

- summary_fn:

  Optional summary function.

- verbose:

  Logical. Default `TRUE`.

## Value

A list of class `sbi_tarp` containing:

- ecp:

  Numeric vector of expected coverage probabilities at each nominal
  level.

- alpha:

  Nominal credibility levels.

- n_sims:

  Number of calibration datasets.

- n_references:

  Number of reference points.

## Details

Implementation of `TARP`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Lemos, P., Coogan, A., Hezaveh, Y. & Perreault-Levasseur, L. (2023).
Sampling-based accuracy testing of posterior estimators for general
inference. *ICML*.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- SBI(simulator, prior, x_obs, n_simulations = 5000)
tarp <- TARP(fit, simulator, prior, n_sims = 200)
plot(fit, type = "tarp", tarp = tarp)
} # }
```
