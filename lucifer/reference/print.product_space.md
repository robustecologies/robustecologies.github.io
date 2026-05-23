# Print method for product_space objects

Prints a concise summary of a product space model selection analysis,
including posterior model probabilities, Bayes factors relative to the
best model, and model indicator mixing diagnostics.

## Usage

``` r
# S3 method for class 'product_space'
print(x, ...)
```

## Arguments

- x:

  An object of class `product_space`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a ProductSpace
model-selection fit produced by
[`ProductSpace`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Carlin, B. P., & Chib, S. (1995). Bayesian model choice via Markov chain
Monte Carlo methods. *Journal of the Royal Statistical Society, Series
B*, 57(3), 473-484.
[doi:10.1111/j.2517-6161.1995.tb02042.x](https://doi.org/10.1111/j.2517-6161.1995.tb02042.x)

## See also

[`ProductSpace`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md),
[`summary.product_space`](https://robustecologies.github.io/lucifer/reference/summary.product_space.md),
[`plot.product_space`](https://robustecologies.github.io/lucifer/reference/plot.product_space.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.product_space
} # }
```
