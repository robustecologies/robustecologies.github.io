# Plot method for `charney_devore_6mode` objects

Plot method for `charney_devore_6mode` objects

## Usage

``` r
# S3 method for class 'charney_devore_6mode'
plot(
  x,
  type = c("timeseries", "phase", "regimes"),
  var = c("x1", "wave"),
  n_points = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `charney_devore_6mode`.

- type:

  Character; the plot to render. One of `"timeseries"`, `"phase"`, or
  `"regimes"`.

- var:

  Character; for `type = "timeseries"`, which observable to plot. One of
  `"x1"` (zonal flow mode) or `"wave"` (square root of wave energy).

- n_points:

  Optional integer; thin the trajectory to this many points before
  rendering. `NULL` plots all points.

- ...:

  Additional graphical parameters passed to the underlying plotting
  backend.

## Value

A ggplot or plotly object (rendered as a side effect, returned
invisibly).

## See also

[`charney_devore_6mode()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_6mode.md),
[`print.charney_devore_6mode()`](https://robustecologies.github.io/tuRbulence/reference/print.charney_devore_6mode.md),
[`summary.charney_devore_6mode()`](https://robustecologies.github.io/tuRbulence/reference/summary.charney_devore_6mode.md),
[`charney_devore_sim()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_sim.md)
