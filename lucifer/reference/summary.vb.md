# Summary method for Variational Bayes results

Provides a structured summary of the Variational Bayes fit including
method, convergence status, ELBO-related diagnostics, and parameter
estimates. When posterior samples are available (via sampling importance
resampling), the summary includes posterior quantiles and effective
sample sizes.

## Usage

``` r
# S3 method for class 'vb'
summary(object, Objective = FALSE, ...)
```

## Arguments

- object:

  An object of class `vb`, as returned by
  [`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

- Objective:

  Logical. If `TRUE`, the deviance trajectory is summarized with
  initial, final, and range statistics. Defaults to `FALSE`.

- ...:

  Additional arguments (currently unused).

## Value

Returns the input object invisibly. The summary is printed as a
side-effect.

## Details

Produces a tabular summary of a Variational Bayes fit produced by
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Blei, D. M., Kucukelbir, A., & McAuliffe, J. D. (2017). Variational
inference: A review for statisticians. *Journal of the American
Statistical Association*, 112(518), 859-877.
[doi:10.1080/01621459.2017.1285773](https://doi.org/10.1080/01621459.2017.1285773)

## See also

[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`print.vb`](https://robustecologies.github.io/lucifer/reference/print.vb.md),
[`plot.vb`](https://robustecologies.github.io/lucifer/reference/plot.vb.md)

## Examples

``` r
if (FALSE) { # \dontrun{
Fit <- VariationalBayes(Model, IV, Data = MyData, Iterations = 1000,
     Method = "Salimans2", Stop.Tolerance = 1e-3)
summary(Fit)
summary(Fit, Objective = TRUE)
} # }
```
