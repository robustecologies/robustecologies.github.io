# Convert external Bayesian fit objects to demonoid

Converts posterior samples from Stan (`stanfit`, `brmsfit`,
`CmdStanMCMC`), JAGS (`mcmc.list`, `runjags`), or the posterior package
(`draws_array`, `draws_matrix`, `draws_df`) into a `demonoid` object
compatible with lucifer's post-processing ecosystem.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `demonoid.default`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `demonoid.mcmc.list`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `demonoid.mcmc`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `demonoid.draws_array`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `demonoid.draws_matrix`.

S3 method: apply [`as()`](https://rdrr.io/r/methods/as.html) to objects
of class `demonoid.draws_df`.

## Usage

``` r
as.demonoid(x, ...)

# Default S3 method
as.demonoid(x, ...)

# S3 method for class 'mcmc.list'
as.demonoid(x, log_lik = NULL, ...)

# S3 method for class 'mcmc'
as.demonoid(x, log_lik = NULL, ...)

# S3 method for class 'draws_array'
as.demonoid(x, log_lik = NULL, ...)

# S3 method for class 'draws_matrix'
as.demonoid(x, log_lik = NULL, ...)

# S3 method for class 'draws_df'
as.demonoid(x, log_lik = NULL, ...)
```

## Arguments

- x:

  See Details.

- ...:

  See Details.

- log_lik:

  An optional \\N \times S\\ matrix of pointwise log-likelihoods (\\N\\
  observations, \\S\\ posterior samples). When provided,
  [`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md),
  [`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
  and
  [`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
  influence analysis become available via the
  [`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md)
  extractor. For `stanfit` objects this can also be a character string
  naming the generated quantity containing pointwise log-likelihoods
  (default `"log_lik"`).

## Value

An object of class `demonoid` with all standard fields populated. The
`$.bridge` element contains provenance metadata.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

The resulting `demonoid` object supports the full diagnostic and
visualization toolkit: [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`Rhat()`](https://robustecologies.github.io/lucifer/reference/Rhat.md),
[`ESS()`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`MCSE()`](https://robustecologies.github.io/lucifer/reference/MCSE.md),
[`caterpillar.plot()`](https://robustecologies.github.io/lucifer/reference/caterpillar.plot.md),
and
[`joint.density.plot()`](https://robustecologies.github.io/lucifer/reference/joint.density.plot.md).
When pointwise log-likelihoods are available, `LOO(log_lik(fit))` and
`WAIC(log_lik(fit))` work directly.

For full compatibility with
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
power-scaling and conflict detection, or with
[`predict()`](https://rdrr.io/r/stats/predict.html), a callable lucifer
Model function and Data list are required. Use
[`lucifer_stan()`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md)
or
[`lucifer_jags()`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md)
to obtain fits with built-in refit capability.

Implementation of `as.demonoid.default`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

Implementation of `as.demonoid.mcmc.list`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

Implementation of `as.demonoid.mcmc`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `as.demonoid.draws_array`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

Implementation of `as.demonoid.draws_matrix`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

Implementation of `as.demonoid.draws_df`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer_stan`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md),
[`lucifer_jags`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md),
[`log_lik`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`to_mcmc_list`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Import a coda mcmc.list
library(coda)
chain1 <- mcmc(matrix(rnorm(2000), 1000, 2,
                      dimnames = list(NULL, c("mu", "sigma"))))
chain2 <- mcmc(matrix(rnorm(2000), 1000, 2,
                      dimnames = list(NULL, c("mu", "sigma"))))
chains <- mcmc.list(chain1, chain2)
fit <- as.demonoid(chains)
summary(fit)
plot(fit)

# Import an rstan fit
library(rstan)
stan_fit <- stan(model_code = "...", data = stan_data)
fit <- as.demonoid(stan_fit)
LOO(log_lik(fit))
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.demonoid.default
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.demonoid.mcmc.list
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.demonoid.mcmc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.demonoid.draws_array
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.demonoid.draws_matrix
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving as.demonoid.draws_df
} # }
```
