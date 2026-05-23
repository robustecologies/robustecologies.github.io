# Leave-one-out cross-validation via Pareto smoothed importance sampling

Computes approximate leave-one-out cross-validation (LOO-CV) using
Pareto smoothed importance sampling (PSIS), as described in Vehtari,
Gelman and Gabry (2017). For each observation, raw importance ratios \\1
/ p(y_i \| \theta^s)\\ are stabilized by fitting a generalized Pareto
distribution to the upper tail and replacing extreme values with
smoothed quantiles. The per-observation Pareto shape parameter \\k\\
diagnoses reliability: \\k \< 0.5\\ indicates reliable estimates, \\0.5
\le k \< 0.7\\ is acceptable, \\0.7 \le k \< 1\\ is problematic, and \\k
\ge 1\\ means the estimate is unreliable.

S3 method: apply `LOO()` to objects of class `default`.

S3 method: apply `LOO()` to objects of class `demonoid`.

S3 method: apply `LOO()` to objects of class `vb`.

S3 method: apply `LOO()` to objects of class `pmc`.

S3 method: apply `LOO()` to objects of class `sde_fit`.

## Usage

``` r
LOO(x, ...)

# Default S3 method
LOO(x, ...)

# S3 method for class 'demonoid'
LOO(x, ...)

# S3 method for class 'vb'
LOO(x, ...)

# S3 method for class 'pmc'
LOO(x, ...)

# S3 method for class 'sde_fit'
LOO(x, ...)
```

## Arguments

- x:

  See Details.

- ...:

  See Details.

## Value

An object of class `lucifer_loo`, which is a list with components:

- elpd_loo:

  The estimated expected log pointwise predictive density.

- p_loo:

  The effective number of parameters.

- looic:

  The LOO information criterion (\\-2 \\
  \mathrm{elpd}\_{\mathrm{loo}}\\).

- se_elpd_loo:

  Standard error of `elpd_loo`.

- se_looic:

  Standard error of `looic`.

- pointwise:

  A data.frame with columns `elpd_loo`, `p_loo`, `looic`, and `pareto_k`
  for each observation.

- pareto_k:

  Numeric vector of Pareto shape parameters.

- diagnostics:

  A list with threshold counts and a reliability assessment.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

LOO-CV avoids the optimism inherent in training-data model assessment by
evaluating each observation under the posterior obtained without that
observation. Direct computation requires refitting the model \\N\\
times, but PSIS provides an efficient approximation from a single
posterior fit. The expected log pointwise predictive density is
\$\$\mathrm{elpd}\_{\mathrm{loo}} = \sum\_{i=1}^{N} \log p(y_i \mid
y\_{-i})\$\$ estimated via importance sampling with Pareto-smoothed
weights. The effective number of parameters is \\p\_{\mathrm{loo}} =
\mathrm{lppd} - \mathrm{elpd}\_{\mathrm{loo}}\\ and the LOO information
criterion is \\\mathrm{looic} = -2 \\ \mathrm{elpd}\_{\mathrm{loo}}\\.

GPD fitting follows Zhang and Stephens (2009) with the weakly
informative prior adjustment of Vehtari et al. (2024). The C++ backend
uses OpenMP parallelization over observations.

Implementation of `LOO.default`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `LOO.demonoid`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `LOO.vb`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `LOO.pmc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `LOO.sde_fit`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). "Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC." *Statistics
and Computing*, 27(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
"Pareto smoothed importance sampling." *Journal of Machine Learning
Research*, 25(72), 1–58.

Zhang, J. and Stephens, M. A. (2009). "A new and efficient estimation
method for the generalized Pareto distribution." *Technometrics*, 51(3),
316–325.

## See also

[`WAIC`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`PSIS`](https://robustecologies.github.io/lucifer/reference/PSIS.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

[`WAIC.default`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`as.demonoid.default`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md).

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`WAIC.demonoid`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.demonoid`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`print.demonoid`](https://robustecologies.github.io/lucifer/reference/print.demonoid.md),
[`summary.demonoid`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.md).

[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`WAIC.vb`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.vb`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.vb`](https://robustecologies.github.io/lucifer/reference/plot.vb.md),
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md),
[`print.vb`](https://robustecologies.github.io/lucifer/reference/print.vb.md),
[`summary.vb`](https://robustecologies.github.io/lucifer/reference/summary.vb.md).

[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`WAIC.pmc`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.pmc`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.pmc`](https://robustecologies.github.io/lucifer/reference/plot.pmc.md),
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md),
[`print.pmc`](https://robustecologies.github.io/lucifer/reference/print.pmc.md),
[`to_mcmc_list.pmc`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md).

[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
[`WAIC.sde_fit`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.sde_fit`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.sde_fit`](https://robustecologies.github.io/lucifer/reference/plot.sde_fit.md),
[`predict.sde_fit`](https://robustecologies.github.io/lucifer/reference/predict.sde_fit.md),
[`print.sde_fit`](https://robustecologies.github.io/lucifer/reference/print.sde_fit.md),
[`summary.sde_fit`](https://robustecologies.github.io/lucifer/reference/summary.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate a pointwise log-likelihood matrix
N <- 50
S <- 1000
log_lik <- matrix(rnorm(N * S, mean = -2, sd = 0.5), nrow = N, ncol = S)

# Compute LOO-PSIS
result <- LOO(log_lik)
print(result)
plot(result)

# Compare with WAIC
WAIC(log_lik)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving LOO.default
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving LOO.demonoid
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving LOO.vb
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving LOO.pmc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving LOO.sde_fit
} # }
```
