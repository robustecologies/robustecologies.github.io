# Plot method for consort objects

Visualizes the diagnostic results from
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## Usage

``` r
# S3 method for class 'consort'
plot(x, type = c("diagnostics", "conditions", "trace", "scorecard"), ...)
```

## Arguments

- x:

  An object of class `consort`.

- type:

  Character. Plot type: `"diagnostics"` (default) shows per-parameter
  ESS and Rhat; `"conditions"` shows a pass/fail bar chart; `"trace"`
  shows the deviance (MCMC) or ELBO (VB) trace from the stored
  diagnostics; `"scorecard"` shows a horizontal gauge of the 0-100
  convergence score with a red-yellow-green gradient.

- ...:

  Currently unused.

## Value

Invisibly returns the ggplot object.

## Details

Produces diagnostic plots of a Consort diagnostic report produced by
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
Yao, Y., Kennedy, L., Gabry, J., Buerkner, P.-C., & Modrak, M. (2020).
Bayesian workflow. arXiv:2011.01808.

## See also

[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`print.consort`](https://robustecologies.github.io/lucifer/reference/print.consort.md),
[`summary.consort`](https://robustecologies.github.io/lucifer/reference/summary.consort.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.consort
} # }
```
