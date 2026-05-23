# Check NUTS diagnostics

Runs all NUTS-specific diagnostics and prints a summary similar to
Stan's `check_hmc_diagnostics`. Checks for divergent transitions, low
E-BFMI (energy Bayesian fraction of missing information), and excessive
tree depth saturation.

## Usage

``` r
check_nuts(fit, ebfmi_threshold = 0.3, max_treedepth_frac = 0.1, warmup = 0L)
```

## Arguments

- fit:

  A `demonoid` object from
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  run with `Algorithm = "NUTS"`.

- ebfmi_threshold:

  Minimum acceptable E-BFMI value. Default 0.3 (Stan's threshold).

- max_treedepth_frac:

  Maximum acceptable fraction of iterations hitting the tree depth
  limit. Default 0.1.

- warmup:

  Integer; number of initial iterations to discard from diagnostic plots
  (not from summary statistics). Default 0.

## Value

An S3 object of class `nuts_check` with `print` and `plot` methods.
Contains diagnostic summaries and per-chain data.

## Details

**Divergent transitions** indicate that the numerical integrator
encountered regions of high curvature where the Hamiltonian energy was
not conserved. Even one divergence suggests the posterior may not be
fully explored. The remedy is typically reparameterization or increasing
`delta` (target acceptance rate) in NUTS Specs.

**E-BFMI** (energy Bayesian fraction of missing information) measures
the efficiency of the momentum resampling step. Low E-BFMI (\< 0.3)
indicates that the sampler has difficulty exploring the energy
distribution, often caused by heavy-tailed posteriors or funnel
geometries.

**Tree depth saturation** occurs when the NUTS tree-building routine
hits its maximum depth repeatedly, preventing the sampler from taking
sufficiently long trajectories. This can reduce effective sample size.
Increasing `Lmax` in Specs may help.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`summary.demonoid`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lucifer(Model, Data, IV, Algorithm = "NUTS",
               Specs = list(A = 1000, delta = 0.8, epsilon = 0.1, Lmax = 10))
diag <- check_nuts(fit)
print(diag)
plot(diag)
} # }
```
