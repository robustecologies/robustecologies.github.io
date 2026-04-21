# Plot method for dde_portrait objects

Renders a DDE portrait as a layered 2D figure (ggplot2) combining frozen
system analysis with actual DDE trajectories. The plot layers correspond
to the components computed by
[`dde_portrait()`](https://robustecologies.github.io/janos/reference/dde_portrait.md):
frozen vector field, nullclines, equilibria (with frozen
classification), manifolds, streamlines, and DDE trajectories.

## Usage

``` r
# S3 method for class 'dde_portrait'
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

  A `dde_portrait` object

- type:

  Character string selecting which components to render. One of
  `"default"` (all computed components), `"field"` (frozen vector field
  only), `"nullclines"`, `"equilibria"`, or `"manifolds"`.

- title:

  Optional plot title. If `NULL`, uses the model name.

- arrow_scale:

  Numeric scaling factor for vector field arrows (default 1).

- arrow_type:

  Arrow head style: `"open"` (default) or `"closed"`.

- feasible:

  Controls shading of infeasible state-space regions. `FALSE` (default)
  disables shading. `TRUE` shades regions where any state variable is
  negative. A named list specifies per-state bounds.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object.

## Details

Equilibrium points are encoded by shape and color following the same
convention as
[`plot.phase_portrait`](https://robustecologies.github.io/janos/reference/plot.phase_portrait.md),
but all labels carry a "(frozen)" suffix to indicate that stability is
assessed from the frozen ODE approximation rather than the full DDE
characteristic spectrum:

- Stable node: filled blue circle

- Unstable node: filled red circle

- Stable focus: filled blue diamond

- Unstable focus: filled red diamond

- Center: open green circle

- Saddle: filled orange triangle

- Saddle-focus: filled orange inverted triangle

When delay spectrum data is available and indicates a stability change
relative to the frozen classification, affected equilibria are marked
with an asterisk in the subtitle and annotated on the plot.

Manifold branches are colored blue (stable) or red (unstable) and
represent the frozen system's separatrices.

DDE trajectories are plotted as solid colored curves, showing the actual
delay-dependent dynamics which may differ substantially from the frozen
flow field.

## Examples

``` r
if (FALSE) { # \dontrun{
pp_dde <- model_spec(
    rhs = list(
        N ~ r * N * (1 - N_delay / K) - a * N * P,
        P ~ b * N * P - d * P
    ),
    delays = list(N_delay = list(state = "N", tau = 2.0)),
    state_names = c("N", "P"),
    parms = list(r = 1.0, K = 10.0, a = 0.1, b = 0.05, d = 0.5),
    init = c(N = 5.0, P = 2.0)
)
dp <- dde_portrait(pp_dde, traj_length = 100)
plot(dp)
plot(dp, type = "field")
plot(dp, type = "nullclines")
plot(dp, type = "equilibria")
} # }
```
