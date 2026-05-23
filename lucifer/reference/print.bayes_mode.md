# Print method for Bayesian mode inference

Displays a concise summary of Bayesian mode inference results in a
bayestestR-inspired tabular format: the posterior distribution over the
number of modes, the probability of multimodality, and mode locations
with credible intervals.

## Usage

``` r
# S3 method for class 'bayes_mode'
print(x, digits = 3, max.modes = 10, ...)
```

## Arguments

- x:

  An object of class `bayes_mode` or `bayes_mode_multi`.

- digits:

  Number of decimal digits. Default 3.

- max.modes:

  Maximum number of mode locations to display. Default 10.

- ...:

  Currently unused.

## Value

Invisibly returns the input object.

## Details

Produces a concise one-screen console report of a Bayesian mode
inference fit produced by
[`BayesMode`](https://robustecologies.github.io/lucifer/reference/BayesMode.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Cross, J. L., Hoogerheide, L., Ulker, P., & van Dijk, H. K. (2024).
Sparse finite mixtures for modal inference. *Economics Letters*, 235,
111551.
[doi:10.1016/j.econlet.2024.111551](https://doi.org/10.1016/j.econlet.2024.111551)

## See also

[`plot.bayes_mode`](https://robustecologies.github.io/lucifer/reference/plot.bayes_mode.md),
[`summary.bayes_mode`](https://robustecologies.github.io/lucifer/reference/summary.bayes_mode.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.bayes_mode
} # }
```
