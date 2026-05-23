# Summary method for SMC objects

Provides an extended summary of an SMC fit, including tempering schedule
details, ESS diagnostics, and marginal likelihood.

## Usage

``` r
# S3 method for class 'smc'
summary(object, ...)
```

## Arguments

- object:

  An object of class `smc`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns a list with summary components.

## Details

Produces a tabular summary of a sequential Monte Carlo fit produced by
[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Del Moral, P., Doucet, A., & Jasra, A. (2006). Sequential Monte Carlo
samplers. *Journal of the Royal Statistical Society, Series B*, 68(3),
411-436.
[doi:10.1111/j.1467-9868.2006.00553.x](https://doi.org/10.1111/j.1467-9868.2006.00553.x)

## See also

[`plot.smc`](https://robustecologies.github.io/lucifer/reference/plot.smc.md),
[`print.smc`](https://robustecologies.github.io/lucifer/reference/print.smc.md),
[`to_mcmc_list.smc`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.smc
} # }
```
