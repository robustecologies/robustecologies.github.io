# Print an object of class `laplace` to the screen

This may be used to print the contents of an object of class `laplace`
to the screen.

## Usage

``` r
# S3 method for class 'laplace'
print(x, ...)
```

## Arguments

- x:

  An object of class `laplace`.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. The side effect is the printed report.

## Details

Produces a concise one-screen console report of a Laplace approximation
fit produced by
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Tierney, L., & Kadane, J. B. (1986). Accurate approximations for
posterior moments and marginal densities. *Journal of the American
Statistical Association*, 81(393), 82-86.
[doi:10.1080/01621459.1986.10478240](https://doi.org/10.1080/01621459.1986.10478240)

## See also

[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## See the LaplaceApproximation function for an example.
} # }
```
