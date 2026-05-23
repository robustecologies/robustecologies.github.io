# Print an object of class `iterquad` to the screen

This may be used to print the contents of an object of class `iterquad`
to the screen.

## Usage

``` r
# S3 method for class 'iterquad'
print(x, ...)
```

## Arguments

- x:

  An object of class `iterquad`.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. The side effect is the printed report.

## Details

Produces a concise one-screen console report of an iterative quadrature
fit produced by
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Naylor, J. C., & Smith, A. F. M. (1982). Applications of a method for
the efficient computation of posterior distributions. *Applied
Statistics*, 31(3), 214-225.
[doi:10.2307/2347995](https://doi.org/10.2307/2347995)

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## See the IterativeQuadrature function for an example.
} # }
```
