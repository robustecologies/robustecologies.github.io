# Convert lucifer fits to coda mcmc.list

Converts a `demonoid`, `vb`, `smc`, or `pmc` object to a
[`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) for use with
coda diagnostics or packages that accept coda objects.

S3 method: apply `to_mcmc_list()` to objects of class `demonoid`.

S3 method: apply `to_mcmc_list()` to objects of class `vb`.

S3 method: apply `to_mcmc_list()` to objects of class `smc`.

S3 method: apply `to_mcmc_list()` to objects of class `pmc`.

## Usage

``` r
to_mcmc_list(x, ...)

# S3 method for class 'demonoid'
to_mcmc_list(x, use_chains = TRUE, ...)

# S3 method for class 'vb'
to_mcmc_list(x, ...)

# S3 method for class 'smc'
to_mcmc_list(x, ...)

# S3 method for class 'pmc'
to_mcmc_list(x, ...)
```

## Arguments

- x:

  A fitted lucifer object.

- ...:

  Additional arguments (currently unused).

- use_chains:

  Logical. If `TRUE` (default) and the fit contains multiple chains,
  each chain is exported as a separate `mcmc` element. If `FALSE`, all
  samples are combined into a single chain.

## Value

An object of class
[`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html).

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `to_mcmc_list`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `to_mcmc_list.demonoid`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

Implementation of `to_mcmc_list.vb`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `to_mcmc_list.smc`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `to_mcmc_list.pmc`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`as.demonoid`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md),
[`to_draws_array`](https://robustecologies.github.io/lucifer/reference/to_draws_array.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LOO.demonoid`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.demonoid`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.demonoid`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`print.demonoid`](https://robustecologies.github.io/lucifer/reference/print.demonoid.md).

[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`LOO.vb`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.vb`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.vb`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.vb`](https://robustecologies.github.io/lucifer/reference/plot.vb.md),
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md),
[`print.vb`](https://robustecologies.github.io/lucifer/reference/print.vb.md).

[`plot.smc`](https://robustecologies.github.io/lucifer/reference/plot.smc.md),
[`print.smc`](https://robustecologies.github.io/lucifer/reference/print.smc.md),
[`summary.smc`](https://robustecologies.github.io/lucifer/reference/summary.smc.md).

[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`LOO.pmc`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.pmc`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.pmc`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.pmc`](https://robustecologies.github.io/lucifer/reference/plot.pmc.md),
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md),
[`print.pmc`](https://robustecologies.github.io/lucifer/reference/print.pmc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lucifer(Model, Data, IV, Algorithm = "NUTS")
ml <- to_mcmc_list(fit)
coda::gelman.diag(ml)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_mcmc_list.demonoid
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_mcmc_list.vb
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_mcmc_list.smc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving to_mcmc_list.pmc
} # }
```
