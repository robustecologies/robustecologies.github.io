# Print method for SMC objects

Prints a concise summary of an SMC fit.

## Usage

``` r
# S3 method for class 'smc'
print(x, ...)
```

## Arguments

- x:

  An object of class `smc`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a sequential Monte Carlo
fit produced by
[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Del Moral, P., Doucet, A., & Jasra, A. (2006). Sequential Monte Carlo
samplers. *Journal of the Royal Statistical Society, Series B*, 68(3),
411-436.
[doi:10.1111/j.1467-9868.2006.00553.x](https://doi.org/10.1111/j.1467-9868.2006.00553.x)

## See also

[`plot.smc`](https://robustecologies.github.io/lucifer/reference/plot.smc.md),
[`summary.smc`](https://robustecologies.github.io/lucifer/reference/summary.smc.md),
[`to_mcmc_list.smc`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.smc
} # }
```
