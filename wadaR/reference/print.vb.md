# Print an object of class `vb` to the screen

This may be used to print the contents of an object of class `vb` to the
screen.

## Usage

``` r
# S3 method for class 'vb'
print(x, ...)
```

## Arguments

- x:

  An object of class `vb`.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. The side effect is the printed report.

## Details

Produces a concise one-screen console report of a Variational Bayes fit
produced by
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Blei, D. M., Kucukelbir, A., & McAuliffe, J. D. (2017). Variational
inference: A review for statisticians. *Journal of the American
Statistical Association*, 112(518), 859-877.
[doi:10.1080/01621459.2017.1285773](https://doi.org/10.1080/01621459.2017.1285773)

## See also

[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## See the VariationalBayes function for an example.
} # }
```
