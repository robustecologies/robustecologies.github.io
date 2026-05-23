# Classifier two-sample test for SBI

Evaluates posterior quality by training a binary classifier to
distinguish samples from the approximate posterior from reference
samples. If the classifier achieves accuracy near 0.5, the two
distributions are indistinguishable (good calibration). Accuracy near
1.0 indicates the distributions differ (poor calibration).

## Usage

``` r
C2ST(
  sbi_fit,
  reference_samples = NULL,
  prior = NULL,
  n_folds = 5L,
  hidden_layers = c(64L, 64L),
  max_epochs = 50L,
  verbose = TRUE
)
```

## Arguments

- sbi_fit:

  An object of class `sbi`.

- reference_samples:

  Matrix of reference samples (e.g., from a known analytical posterior).
  If `NULL`, draws from the prior.

- prior:

  Prior function. Required when `reference_samples` is `NULL`.

- n_folds:

  Integer: number of cross-validation folds. Default 5.

- hidden_layers:

  Integer vector: classifier hidden layers. Default `c(64L, 64L)`.

- max_epochs:

  Integer: maximum training epochs. Default 50.

- verbose:

  Logical. Default `TRUE`.

## Value

A list of class `sbi_c2st` containing:

- accuracy:

  Mean cross-validated classification accuracy.

- accuracies:

  Per-fold accuracies.

- n_folds:

  Number of folds.

## Details

Implementation of `C2ST`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Lopez-Paz, D. & Oquab, M. (2017). Revisiting classifier two-sample
tests. *ICLR*.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- SBI(simulator, prior, x_obs, n_simulations = 5000)
c2st <- C2ST(fit, prior = prior)
print(c2st$accuracy)
} # }
```
