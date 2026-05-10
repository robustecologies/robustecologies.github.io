# Plot method for Bayesian mode inference

Produces publication-quality ggplot2 visualizations of Bayesian mode
inference results using `theme_relab`.

## Usage

``` r
# S3 method for class 'bayes_mode'
plot(
  x,
  type = c("modes", "n_modes", "mixture", "trace", "posterior", "convergence"),
  draws = 100,
  ...
)
```

## Arguments

- x:

  An object of class `bayes_mode` or `bayes_mode_multi`.

- type:

  Character string specifying the plot type: `"modes"` (default) shows
  data density with posterior mode locations as colored vertical bands;
  `"n_modes"` shows the posterior distribution over the number of modes;
  `"mixture"` shows the fitted mixture density with individual
  components; `"trace"` shows MCMC trace plots for alpha, effective K,
  and number of modes; `"posterior"` shows violin plots of mode location
  posteriors; `"convergence"` shows MCMC convergence diagnostics.

- draws:

  Number of mixture draws to overlay in `"mixture"` plots. Default 100.

- ...:

  Currently unused.

## Value

A `ggplot` object (invisibly).

## Details

Produces diagnostic plots of a Bayesian mode inference fit produced by
[`BayesMode`](https://robustecologies.github.io/lucifer/reference/BayesMode.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Cross, J. L., Hoogerheide, L., Ulker, P., & van Dijk, H. K. (2024).
Sparse finite mixtures for modal inference. *Economics Letters*, 235,
111551.
[doi:10.1016/j.econlet.2024.111551](https://doi.org/10.1016/j.econlet.2024.111551)

## See also

[`print.bayes_mode`](https://robustecologies.github.io/lucifer/reference/print.bayes_mode.md),
[`summary.bayes_mode`](https://robustecologies.github.io/lucifer/reference/summary.bayes_mode.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.bayes_mode
} # }
```
