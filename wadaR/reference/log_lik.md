# Extract pointwise log-likelihoods from a fitted model

Returns the \\N \times S\\ matrix of pointwise log-likelihoods stored in
a fitted model object. This matrix is the required input for
[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md) and
[`WAIC`](https://robustecologies.github.io/lucifer/reference/WAIC.md).

S3 method: apply `log_lik()` to objects of class `demonoid`.

S3 method: apply `log_lik()` to objects of class `vb`.

S3 method: apply `log_lik()` to objects of class `pmc`.

S3 method: apply `log_lik()` to objects of class `sde_fit`.

## Usage

``` r
log_lik(x, ...)

# S3 method for class 'demonoid'
log_lik(x, ...)

# S3 method for class 'vb'
log_lik(x, ...)

# S3 method for class 'pmc'
log_lik(x, ...)

# S3 method for class 'sde_fit'
log_lik(x, ...)
```

## Arguments

- x:

  A fitted model object (class `demonoid`, `vb`, `smc`, or `pmc`).

- ...:

  Additional arguments (currently unused).

## Value

An \\N \times S\\ matrix of pointwise log-likelihoods.

See Details.

See Details.

See Details.

See Details.

## Details

For `demonoid` objects created via
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
or
[`lucifer_stan()`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md)/[`lucifer_jags()`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md),
the log-likelihood matrix is stored in the `$log_lik` field. For native
lucifer fits where the Monitor matrix contains pointwise log-likelihoods
(one column per observation), `log_lik()` returns `t(Monitor)`.

Implementation of `log_lik.demonoid`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `log_lik.vb`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `log_lik.pmc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `log_lik.sde_fit`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`as.demonoid`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LOO.demonoid`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.demonoid`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`print.demonoid`](https://robustecologies.github.io/lucifer/reference/print.demonoid.md),
[`summary.demonoid`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.md).

[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`LOO.vb`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.vb`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`plot.vb`](https://robustecologies.github.io/lucifer/reference/plot.vb.md),
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md),
[`print.vb`](https://robustecologies.github.io/lucifer/reference/print.vb.md),
[`summary.vb`](https://robustecologies.github.io/lucifer/reference/summary.vb.md).

[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`LOO.pmc`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.pmc`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`plot.pmc`](https://robustecologies.github.io/lucifer/reference/plot.pmc.md),
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md),
[`print.pmc`](https://robustecologies.github.io/lucifer/reference/print.pmc.md),
[`to_mcmc_list.pmc`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md).

[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
[`LOO.sde_fit`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.sde_fit`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`plot.sde_fit`](https://robustecologies.github.io/lucifer/reference/plot.sde_fit.md),
[`predict.sde_fit`](https://robustecologies.github.io/lucifer/reference/predict.sde_fit.md),
[`print.sde_fit`](https://robustecologies.github.io/lucifer/reference/print.sde_fit.md),
[`summary.sde_fit`](https://robustecologies.github.io/lucifer/reference/summary.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# After importing a Stan fit with log_lik
fit <- as.demonoid(stan_fit)
ll <- log_lik(fit)
LOO(ll)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving log_lik.demonoid
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving log_lik.vb
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving log_lik.pmc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving log_lik.sde_fit
} # }
```
