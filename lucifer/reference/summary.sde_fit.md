# Summary method for SDE fit objects

Produces an extended summary of an SDE model fit, including parameter
posterior distributions, MCMC diagnostics, and model fit metrics.

## Usage

``` r
# S3 method for class 'sde_fit'
summary(object, ...)
```

## Arguments

- object:

  An object of class `sde_fit`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns the summary information as a list.

## Details

Produces a tabular summary of a stochastic differential equation fit
produced by
[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Iacus, S. M. (2008). *Simulation and Inference for Stochastic
Differential Equations*. Springer.
[doi:10.1007/978-0-387-75839-8](https://doi.org/10.1007/978-0-387-75839-8)

## See also

[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
[`LOO.sde_fit`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.sde_fit`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.sde_fit`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.sde_fit`](https://robustecologies.github.io/lucifer/reference/plot.sde_fit.md),
[`predict.sde_fit`](https://robustecologies.github.io/lucifer/reference/predict.sde_fit.md),
[`print.sde_fit`](https://robustecologies.github.io/lucifer/reference/print.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.sde_fit
} # }
```
