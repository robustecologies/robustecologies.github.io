# Plot method for mode estimates

Visualizes mode estimation results using ggplot2 with `theme_relab`.
Three plot types are available: density overlay, histogram with markers,
and cross-method comparison.

## Usage

``` r
# S3 method for class 'mode_estimate'
plot(x, type = c("density", "histogram", "comparison"), bins = 30, ...)
```

## Arguments

- x:

  An object of class `mode_estimate`.

- type:

  Character string: `"density"` (default) shows a kernel density curve
  with vertical lines at detected modes; `"histogram"` shows a histogram
  with mode markers; `"comparison"` shows a faceted panel comparing
  modes from all four estimation methods.

- bins:

  Number of histogram bins. Default 30.

- ...:

  Currently unused.

## Value

A `ggplot` object (invisibly).

## Details

Produces diagnostic plots of a mode estimate produced by
[`Mode`](https://robustecologies.github.io/lucifer/reference/Mode.md).
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

[`Mode`](https://robustecologies.github.io/lucifer/reference/Mode.md),
[`as.double.mode_estimate`](https://robustecologies.github.io/lucifer/reference/as.double.mode_estimate.md),
[`print.mode_estimate`](https://robustecologies.github.io/lucifer/reference/print.mode_estimate.md),
[`summary.mode_estimate`](https://robustecologies.github.io/lucifer/reference/summary.mode_estimate.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.mode_estimate
} # }
```
