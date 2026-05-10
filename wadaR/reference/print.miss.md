# Print an object of class `miss` to the screen

This may be used to print the contents of an object of class `miss` to
the screen.

## Usage

``` r
# S3 method for class 'miss'
print(x, ...)
```

## Arguments

- x:

  An object of class `miss`.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. The side effect is the printed report.

## Details

Produces a concise one-screen console report of a multiple-imputation
fit produced by
[`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Little, R. J. A., & Rubin, D. B. (2019). *Statistical Analysis with
Missing Data* (3rd ed.). Wiley. ISBN 9780470526798.

## See also

[`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## See the MISS function for an example.
} # }
```
