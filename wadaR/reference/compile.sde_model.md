# Compile SDE expressions to native C++

Translates drift and diffusion expressions from R to C++, compiles at
runtime via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html), and
extracts native function pointers for zero-overhead dispatch in the
particle filter and EKF engines. The compiled code is cached based on a
hash of the expressions, so repeated calls with the same model skip
compilation.

## Usage

``` r
# S3 method for class 'sde_model'
compile(spec, Data = NULL, cache = TRUE, verbose = TRUE, ...)
```

## Arguments

- spec:

  An object of class `sde_model` whose `drift` and `diffusion` were
  specified as character strings, expressions, or formulas (not R
  functions).

- Data:

  Ignored. Accepted for signature compatibility with the
  `compile.model_spec` method.

- cache:

  Ignored. Accepted for signature compatibility; SDE compilation always
  caches based on a content hash of the expressions.

- verbose:

  Logical; if `TRUE`, print compilation messages.

- ...:

  Additional arguments (currently unused).

## Value

The input `sde_model` object with native function pointers stored in
`object$Data$.native_drift_ptr` and `object$Data$.native_diffusion_ptr`,
and `object$Data$.sde_compiled = TRUE`.

## Details

Compile an sde_model's drift/diffusion expressions to C++

Implementation of `compile.sde_model`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md),
[`plot.sde_model`](https://robustecologies.github.io/lucifer/reference/plot.sde_model.md),
[`print.sde_model`](https://robustecologies.github.io/lucifer/reference/print.sde_model.md),
[`simulate.sde_model`](https://robustecologies.github.io/lucifer/reference/simulate.sde_model.md),
[`summary.sde_model`](https://robustecologies.github.io/lucifer/reference/summary.sde_model.md).

## Examples

``` r
if (FALSE) { # \dontrun{
sde <- SDE(
    drift = "theta[1] * (theta[2] - x[1])",
    diffusion = "theta[3]",
    data = y, times = times, x0 = y[1],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[4],
    prior = function(theta) sum(dnorm(theta, 0, 10, log = TRUE)),
    parm.names = c("kappa", "mu", "sigma", "obs_sd"),
    compile = TRUE
)
} # }
```
