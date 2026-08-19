# Plot a multi-start ensemble

Views of the start-level outcomes of a multi-start run.

## Usage

``` r
# S3 method for class 'sym_ensemble'
plot(x, type = c("starts", "basins", "spread", "histogram", "dashboard"), ...)
```

## Arguments

- x:

  A `sym_ensemble` object.

- type:

  View: `"starts"` (incumbent value per start, ordered), `"basins"`
  (start points colored by their attained value, dimension 2 only),
  `"spread"` (attained value against the distance from the start to the
  ensemble incumbent, any dimension), `"histogram"` (distribution of the
  attained values), or `"dashboard"` (a four-panel composite).

- ...:

  Ignored.

## Value

The ggplot object, invisibly; for `type = "dashboard"` the assembled
grid graphical object.

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`print.sym_ensemble()`](https://robustecologies.github.io/symplectoR/reference/print.sym_ensemble.md),
[`summary.sym_ensemble()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_ensemble.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ens <- sym_optim(sym_objective(sym_benchmark("rosenbrock", d = 2)),
                 x0 = c(0, 0), method = "rgd", n_starts = 16, seed = 1)
plot(ens, type = "starts")
plot(ens, type = "basins")
plot(ens, type = "dashboard")
} # }
```
