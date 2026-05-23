# Print method for SDE fit objects

Prints a concise summary of an SDE model fit.

## Usage

``` r
# S3 method for class 'sde_fit'
print(x, ...)
```

## Arguments

- x:

  An object of class `sde_fit`.

- ...:

  Additional arguments passed to `print.demonoid`.

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a stochastic
differential equation fit produced by
[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

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
[`summary.sde_fit`](https://robustecologies.github.io/lucifer/reference/summary.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.sde_fit
} # }
```
