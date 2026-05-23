# Summary method for a correlation_dimension object

Reports the Grassberger-Procaccia correlation dimension \\D_2\\ with its
95% confidence interval from the scaling-window linear fit, the
coefficient of determination \\R^2\\ of the fit, the number of
\\\varepsilon\\ radii inside the scaling window, the Theiler window, the
sample size and the fitted scaling range in log units.

## Usage

``` r
# S3 method for class 'correlation_dimension'
summary(object, ...)
```

## Arguments

- object:

  A `correlation_dimension` object.

- ...:

  Unused.

## Value

A list of class `summary.correlation_dimension`, invisibly.

## Details

The confidence interval is computed via
`stats::confint(x$fit, level = 0.95)` on the slope coefficient of the
linear regression of \\\log C(\varepsilon)\\ against \\\log
\varepsilon\\ restricted to the scaling window. A narrow interval with
\\R^2 \approx 1\\ indicates a well-resolved power-law scaling; a wide
interval or low \\R^2\\ suggests too small a sample, too narrow a
scaling window, or residual non-stationarity.

## References

Grassberger, P., & Procaccia, I. (1983). Measuring the strangeness of
strange attractors. *Physica D*, 9(1-2), 189-208.
[doi:10.1016/0167-2789(83)90298-1](https://doi.org/10.1016/0167-2789%2883%2990298-1)

Theiler, J. (1986). Spurious dimension from correlation algorithms
applied to limited time-series data. *Physical Review A*, 34(3),
2427-2432.
[doi:10.1103/PhysRevA.34.2427](https://doi.org/10.1103/PhysRevA.34.2427)

## See also

[`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)
. constructor;
[`print.correlation_dimension()`](https://robustecologies.github.io/janos/reference/print.correlation_dimension.md)
. compact header;
[`plot.correlation_dimension()`](https://robustecologies.github.io/janos/reference/plot.correlation_dimension.md)
. log-log scaling plot.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50)
summary(correlation_dimension(run))
} # }
```
