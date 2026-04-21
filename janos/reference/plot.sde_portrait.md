# Plot method for sde_portrait objects

Renders a stochastic phase portrait as a layered 2D figure (ggplot2) or
an interactive 3D scene (plotly). The plot layers combine the
deterministic skeleton structure with stochastic elements: diffusion
intensity heatmap, sample paths, and confidence ellipses.

## Usage

``` r
# S3 method for class 'sde_portrait'
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

  An `sde_portrait` object

- type:

  Character string selecting which components to render. One of
  `"default"`, `"drift"`, `"diffusion"`, `"skeleton"`, or
  `"stochastic"`.

- title:

  Optional plot title. If `NULL`, uses the model name.

- arrow_scale:

  Numeric scaling factor for drift field arrows (default 1).

- arrow_type:

  Arrow head style: `"open"` (default) or `"closed"`.

- feasible:

  Controls shading of infeasible regions (2D only). `FALSE` (default)
  disables shading. `TRUE` shades regions where any state variable is
  negative. A named list specifies per-state bounds. Ignored for 3D
  portraits.

- ...:

  Additional arguments passed to underlying plot functions.

## Value

For 2D: a ggplot2 object. For 3D: a plotly object.

## Details

The `type` argument selects which layers are rendered:

- `"default"`: all computed components.

- `"drift"`: drift field only.

- `"diffusion"`: diffusion intensity heatmap only.

- `"skeleton"`: deterministic skeleton only (drift, nullclines,
  equilibria, manifolds). No noise-related layers.

- `"stochastic"`: stochastic layers only (sample paths, confidence
  ellipses, diffusion heatmap), plus equilibria for context.

Equilibrium points are encoded by shape and color following the same
convention as `plot.phase_portrait`:

- Stable node: filled blue circle

- Unstable node: filled red circle

- Stable focus: filled blue diamond

- Unstable focus: filled red diamond

- Center: open green circle

- Saddle: filled orange triangle

- Saddle-focus: filled orange inverted triangle

## Examples

``` r
if (FALSE) { # \dontrun{
ou <- model_spec(
    rhs = list(x ~ -a * x, y ~ -b * y),
    diffusion = list(x ~ sigma, y ~ sigma),
    state_names = c("x", "y"),
    parms = list(a = 1, b = 2, sigma = 0.3),
    init = c(x = 1, y = 1)
)
set.seed(42)
sp <- sde_portrait(ou, xlim = c(-3, 3), ylim = c(-3, 3))

plot(sp)
plot(sp, type = "drift")
plot(sp, type = "diffusion")
plot(sp, type = "skeleton")
plot(sp, type = "stochastic")
} # }
```
