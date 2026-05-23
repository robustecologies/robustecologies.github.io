# Plot method for a correlation_dimension object

Renders the log-log scaling plot of the Grassberger-Procaccia
correlation sum \\C(\varepsilon)\\ against the radius \\\varepsilon\\.
Points inside the scaling window used for the \\D_2\\ fit are coloured
differently from points outside, and the fitted line of slope \\D_2\\ is
overlaid for visual inspection of the power-law regime.

## Usage

``` r
# S3 method for class 'correlation_dimension'
plot(x, title = NULL, ...)
```

## Arguments

- x:

  A `correlation_dimension` object.

- title:

  Optional plot title (overrides the default).

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot object.

## See also

[`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)
. constructor;
[`print.correlation_dimension()`](https://robustecologies.github.io/janos/reference/print.correlation_dimension.md)
. compact header;
[`summary.correlation_dimension()`](https://robustecologies.github.io/janos/reference/summary.correlation_dimension.md)
. fit diagnostics.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50)
plot(correlation_dimension(run))
} # }
```
