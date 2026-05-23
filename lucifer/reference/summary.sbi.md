# Summary method for SBI objects

Provides an extended summary of a simulation-based inference fit,
including network architecture, training diagnostics, and posterior
statistics. For NLE and NRE methods, MCMC convergence diagnostics are
included.

## Usage

``` r
# S3 method for class 'sbi'
summary(object, ...)
```

## Arguments

- object:

  An object of class `sbi`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns a list with summary components.

## Details

Produces a tabular summary of a simulation-based inference fit produced
by [`SBI`](https://robustecologies.github.io/lucifer/reference/SBI.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Cranmer, K., Brehmer, J., & Louppe, G. (2020). The frontier of
simulation-based inference. *PNAS*, 117(48), 30055-30062.
[doi:10.1073/pnas.1912789117](https://doi.org/10.1073/pnas.1912789117)

## See also

[`plot.sbi`](https://robustecologies.github.io/lucifer/reference/plot.sbi.md),
[`print.sbi`](https://robustecologies.github.io/lucifer/reference/print.sbi.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.sbi
} # }
```
