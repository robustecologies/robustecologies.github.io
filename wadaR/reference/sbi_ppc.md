# Posterior predictive checks for SBI

Generates replicated data from the posterior predictive distribution by
sampling parameters from the SBI posterior and simulating new datasets.
Useful for assessing model fit.

## Usage

``` r
sbi_ppc(sbi_fit, simulator, n_sims = 100L, summary_fn = NULL)
```

## Arguments

- sbi_fit:

  An object of class `sbi`.

- simulator:

  The simulator function.

- n_sims:

  Integer: number of posterior predictive replications. Default 100.

- summary_fn:

  Optional summary function.

## Value

A list of length `n_sims`, each element containing a simulated data
vector (or summary statistics).

## Details

Implementation of `sbi_ppc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- SBI(simulator, prior, x_obs, n_simulations = 5000)
ppc <- sbi_ppc(fit, simulator)
plot(fit, type = "ppc", ppc_data = ppc)
} # }
```
