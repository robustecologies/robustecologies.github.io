# Print an object of class `bayesquad`

Prints a concise summary of a `bayesquad` object produced by
[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md).

## Usage

``` r
# S3 method for class 'bayesquad'
print(x, ...)
```

## Arguments

- x:

  An object of class `bayesquad`.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. The side effect is the printed report.

## Details

Produces a concise one-screen console report of a Bayesian quadrature
fit produced by
[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

O'Hagan, A. (1991). Bayes-Hermite quadrature. *Journal of Statistical
Planning and Inference*, 29(3), 245-260.
[doi:10.1016/0378-3758(91)90002-V](https://doi.org/10.1016/0378-3758%2891%2990002-V)

## See also

[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## See the BayesianQuadrature function for an example.
} # }
```
