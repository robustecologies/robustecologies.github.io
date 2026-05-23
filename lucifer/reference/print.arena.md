# Print an arena object

Displays a compact metrics table from the universal benchmarking
comparison, highlighting the best method and pass/fail convergence
status.

## Usage

``` r
# S3 method for class 'arena'
print(x, ...)
```

## Arguments

- x:

  An object of class `arena`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of an Arena competition
result produced by
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman &
Hall/CRC. ISBN 9781439840955.

## See also

[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`summary.arena`](https://robustecologies.github.io/lucifer/reference/summary.arena.md),
[`plot.arena`](https://robustecologies.github.io/lucifer/reference/plot.arena.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.arena
} # }
```
