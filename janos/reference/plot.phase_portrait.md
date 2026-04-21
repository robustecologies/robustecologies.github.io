# Plot method for phase_portrait objects

Renders a phase portrait as a layered 2D figure (ggplot2) or an
interactive 3D scene (plotly). The plot layers correspond to the
components computed by
[`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md):
vector field, nullclines, equilibria, streamlines, manifolds, and
trajectories.

## Usage

``` r
# S3 method for class 'phase_portrait'
plot(
  x,
  type = "default",
  title = NULL,
  arrow_scale = 1,
  arrow_type = c("open", "closed"),
  feasible = FALSE,
  ...
)
```

## Arguments

- x:

  A `phase_portrait` object

- type:

  Character string selecting which components to render. One of
  `"default"` (all computed components), `"field"` (vector field only),
  `"nullclines"`, `"streamlines"`, `"equilibria"`, or `"manifolds"`.

- title:

  Optional plot title. If `NULL`, uses the model name.

- arrow_scale:

  Numeric scaling factor for vector field arrows (default 1).

- arrow_type:

  Arrow head style: `"open"` (default) draws a V-shaped head that does
  not obscure the shaft, while `"closed"` draws a filled triangular
  head.

- feasible:

  Controls shading of infeasible state-space regions (2D only). `FALSE`
  (default) disables shading. `TRUE` shades regions where any state
  variable is negative. A named list (e.g.
  `list(x = c(0, 100), y = c(0, 50))`) shades regions outside the
  specified per-state bounds; states not listed are unconstrained.
  Ignored for 3D portraits.

- ...:

  Additional arguments passed to underlying plot functions.

## Value

For 2D: a ggplot2 object. For 3D: a plotly object.

## Details

Equilibrium points are encoded by both shape and color to remain
distinguishable for colorblind readers:

- Stable node: filled blue circle

- Unstable node: filled red circle

- Stable focus: filled blue diamond

- Unstable focus: filled red diamond

- Center: open green circle

- Saddle: filled orange triangle

- Saddle-focus: filled orange inverted triangle

Manifold branches are colored blue (stable) or red (unstable).

## Examples

``` r
if (FALSE) { # \dontrun{
lv <- model_spec(
    rhs = list(x ~ a * x - b * x * y, y ~ d * x * y - c * y),
    state_names = c("x", "y"),
    parms = list(a = 1.0, b = 0.1, d = 0.075, c = 1.5),
    init = c(x = 20, y = 10)
)
pp <- phase_portrait(lv, xlim = c(0, 40), ylim = c(0, 25),
                     streamlines = TRUE)

plot(pp)
plot(pp, type = "field")
plot(pp, type = "nullclines")
plot(pp, type = "equilibria")
} # }
```
