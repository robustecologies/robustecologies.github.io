# Summarize an arena object

Provides the full metrics table, pairwise comparison matrices, Pareto
frontier, per-parameter mean/SD comparison, and method-specific
diagnostics.

## Usage

``` r
# S3 method for class 'arena'
summary(object, Parms = NULL, ...)
```

## Arguments

- object:

  An object of class `arena`.

- Parms:

  Integer vector or character vector of parameter indices or names to
  include in the per-parameter comparison. Default `NULL` shows all (up
  to 20).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `object`.

## Details

Produces a tabular summary of an Arena competition result produced by
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman &
Hall/CRC. ISBN 9781439840955.

## See also

[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`print.arena`](https://robustecologies.github.io/lucifer/reference/print.arena.md),
[`plot.arena`](https://robustecologies.github.io/lucifer/reference/plot.arena.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.arena
} # }
```
