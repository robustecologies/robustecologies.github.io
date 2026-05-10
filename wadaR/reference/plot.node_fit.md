# Plot method for node_fit objects

Produces diagnostic and analysis plots for a neural ODE fit. Eight plot
types are available, selected via the `type` argument.

## Usage

``` r
# S3 method for class 'node_fit'
plot(x, type = "interpolation", col = NULL, ...)
```

## Arguments

- x:

  An object of class `node_fit`.

- type:

  Character. One of `"interpolation"` (default), `"dynamics"`,
  `"effects"`, `"contributions"`, `"network"`, `"phase"`, `"cv"`, or
  `"forecast"`.

- col:

  Optional character vector of colors.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns `NULL`.

## Details

Produces diagnostic plots of a neural ODE fit produced by
[`NODE`](https://robustecologies.github.io/lucifer/reference/NODE.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Bonnaffe, W., Sheldon, B. C., & Coulson, T. (2023). Neural ordinary
differential equations for ecological and evolutionary time series
analysis. *Methods in Ecology and Evolution*, 12(7), 1301-1315.
[doi:10.1111/2041-210X.13606](https://doi.org/10.1111/2041-210X.13606)

## See also

[`NODE`](https://robustecologies.github.io/lucifer/reference/NODE.md),
[`predict.node_fit`](https://robustecologies.github.io/lucifer/reference/NODE_predict.md),
[`print.node_fit`](https://robustecologies.github.io/lucifer/reference/print.node_fit.md),
[`summary.node_fit`](https://robustecologies.github.io/lucifer/reference/summary.node_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.node_fit
} # }
```
