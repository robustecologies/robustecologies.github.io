# Plot an objective specification

Surface, contour and axis-profile views of an objective over its box (or
a default window).

## Usage

``` r
# S3 method for class 'sym_objective'
plot(x, type = c("contour", "surface", "profiles", "dashboard"), n = 101, ...)
```

## Arguments

- x:

  A `sym_objective` or `sym_benchmark` object.

- type:

  View: `"contour"` (filled contours, dimension 2; profile curve in
  dimension 1), `"surface"` (raster, dimension 2 only), `"profiles"`
  (the objective along each coordinate axis through the reference point,
  any dimension), or `"dashboard"` (a composite of the views available
  at this dimension).

- n:

  Grid resolution per axis.

- ...:

  Ignored.

## Value

The ggplot object, invisibly; for `type = "dashboard"` the assembled
grid graphical object.

## See also

[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md),
[`print.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/print.sym_objective.md),
[`summary.sym_objective()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_objective.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot(sym_objective(sym_benchmark("rosenbrock", d = 2), lower = -2, upper = 2),
     type = "contour")
plot(sym_objective(sym_benchmark("rosenbrock", d = 2), lower = -2, upper = 2),
     type = "dashboard")
} # }
```
