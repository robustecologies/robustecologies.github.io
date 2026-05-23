# Plot method for crucible objects

Visualizes the automated inference pipeline results.

## Usage

``` r
# S3 method for class 'crucible'
plot(
  x,
  type = c("dashboard", "prescribe", "arena", "convergence", "accuracy", "diagnostics",
    "refinement"),
  ...
)
```

## Arguments

- x:

  An object of class `crucible`.

- type:

  Character. Plot type: `"dashboard"` (default) shows a 2x2 overview;
  `"prescribe"` delegates to `plot(x$Prescription)`; `"arena"` delegates
  to `plot(x$Arena)`; `"convergence"` shows ESS trajectory across
  rounds; `"accuracy"` delegates to `plot(x$Arena, type = "accuracy")`;
  `"diagnostics"` shows the best method's Consort diagnostics;
  `"refinement"` shows a condition-by-round heatmap per method,
  revealing which diagnostics drive refinement iterations.

- ...:

  Passed to downstream plot methods.

## Value

Invisibly returns the plot object.

## Details

Produces diagnostic plots of a Crucible pipeline fit produced by
[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md).
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

[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md),
[`print.crucible`](https://robustecologies.github.io/lucifer/reference/print.crucible.md),
[`summary.crucible`](https://robustecologies.github.io/lucifer/reference/summary.crucible.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.crucible
} # }
```
