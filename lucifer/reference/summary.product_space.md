# Summary method for product_space objects

Provides an extended summary of a product space model selection
analysis, including the full Bayes factor matrix, transition matrix,
per-model posterior summaries, convergence diagnostics, and bisection
calibration history.

## Usage

``` r
# S3 method for class 'product_space'
summary(object, ...)
```

## Arguments

- object:

  An object of class `product_space`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns a list with summary components.

## Details

Produces a tabular summary of a ProductSpace model-selection fit
produced by
[`ProductSpace`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

Carlin, B. P., & Chib, S. (1995). Bayesian model choice via Markov chain
Monte Carlo methods. *Journal of the Royal Statistical Society, Series
B*, 57(3), 473-484.
[doi:10.1111/j.2517-6161.1995.tb02042.x](https://doi.org/10.1111/j.2517-6161.1995.tb02042.x)

## See also

[`ProductSpace`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md),
[`print.product_space`](https://robustecologies.github.io/lucifer/reference/print.product_space.md),
[`plot.product_space`](https://robustecologies.github.io/lucifer/reference/plot.product_space.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.product_space
} # }
```
