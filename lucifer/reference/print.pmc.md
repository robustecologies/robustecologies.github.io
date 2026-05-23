# Print an object of class `pmc` to the screen

This may be used to print the contents of an object of class `pmc` to
the screen.

## Usage

``` r
# S3 method for class 'pmc'
print(x, ...)
```

## Arguments

- x:

  An object of class `pmc`.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. The side effect is the printed report.

## Details

Produces a concise one-screen console report of a Population Monte Carlo
fit produced by
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Cappe, O., Guillin, A., Marin, J. M., & Robert, C. P. (2004). Population
Monte Carlo. *Journal of Computational and Graphical Statistics*, 13(4),
907-929.
[doi:10.1198/106186004X12803](https://doi.org/10.1198/106186004X12803)

## See also

[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## See the PMC function for an example.
} # }
```
