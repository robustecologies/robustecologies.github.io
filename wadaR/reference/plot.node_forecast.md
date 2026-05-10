# Plot method for node_forecast objects

Produces publication-quality forecast visualizations. The default type
overlays observed data, fitted trajectory, and forecast trajectory with
a vertical separator at the forecast origin.

## Usage

``` r
# S3 method for class 'node_forecast'
plot(x, type = "trajectory", col = NULL, ...)
```

## Arguments

- x:

  An object of class `node_forecast`.

- type:

  Character. One of `"trajectory"` (default) or `"phase"`.

- col:

  Optional character vector of colors overriding the default palette.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns `NULL`.

## Details

- `"trajectory"`:

  Observed data (grey points), fitted values (solid colored line), and
  forecast (dashed colored line) for each variable in a faceted layout.
  A vertical dashed line marks the forecast origin. An annotated region
  labels the forecast zone.

- `"phase"`:

  Phase space trajectories for all variable pairs. Observed trajectory
  in grey, fitted in solid color, forecast in dashed color. A marker
  indicates the forecast start point.

## See also

[`print.node_forecast`](https://robustecologies.github.io/lucifer/reference/print.node_forecast.md).

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- NODE_predict(fit, horizon = 5)
plot(pred)
plot(pred, type = "phase")
} # }
```
