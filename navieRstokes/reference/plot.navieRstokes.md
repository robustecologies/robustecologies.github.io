# Plot method for navieRstokes objects

S3 plot method for objects of class 'navieRstokes' produced by
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md).
Returns a ggplot2 object (or a patchwork composition for
`plot_type = "all"`) that can be composed, faceted or saved like any
other ggplot2 plot.

## Usage

``` r
# S3 method for class 'navieRstokes'
plot(
  x,
  time_index = NULL,
  plot_type = c("vorticity", "pressure", "speed", "velocity", "all"),
  subsample = 4L,
  scale = 0.04,
  n_contours = 12L,
  ...
)
```

## Arguments

- x:

  Object of class 'navieRstokes' (output from
  [`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)).

- time_index:

  Integer. Time-step index to render. Defaults to the last saved
  snapshot.

- plot_type:

  One of `"vorticity"`, `"pressure"`, `"speed"`, `"velocity"`, `"all"`.
  Selects the scalar field being coloured. `"velocity"` overlays a
  quiver of velocity vectors on the speed field. `"all"` returns a
  three-panel patchwork composition (vorticity \| speed \| pressure).

- subsample:

  Integer. Stride used when subsampling the velocity grid for the quiver
  overlay (`plot_type = "velocity"`). One arrow is drawn per
  `subsample`-by-`subsample` cell block.

- scale:

  Numeric. Length scaling for the velocity arrows.

- n_contours:

  Integer. Number of contour lines overlaid on the vorticity and
  pressure fields. Set to `0L` to suppress contours. Ignored for
  `"speed"` and `"velocity"`.

- ...:

  Currently unused.

## Value

A ggplot2 object, or a patchwork composition when `plot_type = "all"`
(invisibly when called for its side effect).

## Details

Field colouring uses
[`ggplot2::scale_fill_gradient2`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)
(diverging blue-white-red, midpoint 0) for vorticity,
[`ggplot2::scale_fill_distiller`](https://ggplot2.tidyverse.org/reference/scale_brewer.html)
with palette `"YlOrRd"` for pressure, and
`viridis::scale_fill_viridis_c` for speed and velocity. Grid cells are
rendered at their exact discrete values (`interpolate = FALSE`);
iso-contour lines (`geom_contour`) are overlaid on vorticity and
pressure plots for `n_contours > 0`. Greek labels in the title are
rendered through
[`ggplot2::labs`](https://ggplot2.tidyverse.org/reference/labs.html);
the package recommends
[`ragg::agg_png()`](https://ragg.r-lib.org/reference/agg_png.html) for
saving figures so that Greek glyphs render cleanly on common locales.

## See also

[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md),
[`flow`](https://robustecologies.github.io/navieRstokes/reference/flow.md),
[`compute_vorticity`](https://robustecologies.github.io/navieRstokes/reference/compute_vorticity.md),
[`print.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/print.navieRstokes.md),
[`summary.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/summary.navieRstokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1, ly = 1,
  dt = 0.001, nt = 1000, nu = 0.01,
  initial_condition = function(x, y) {
    rotating_ic(x, y, omega = 2 * pi, x0 = 0.5, y0 = 0.5)
  },
  boundary_condition = list(
    type = "dirichlet",
    values = list(
      u_left = 0, u_right = 0, u_top = 0, u_bottom = 0,
      v_left = 0, v_right = 0, v_top = 0, v_bottom = 0
    )
  ),
  save_interval = 10
)

plot(result, plot_type = "vorticity")
plot(result, plot_type = "velocity")
plot(result, plot_type = "pressure")
plot(result, plot_type = "all")
} # }
```
