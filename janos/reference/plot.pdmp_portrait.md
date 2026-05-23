# Plot method for pdmp_portrait objects

Renders a PDMP portrait as a layered 2D figure using ggplot2. The
default plot shows a faceted view with one panel per regime, displaying
each regime's vector field, nullclines, and equilibria. Sample switching
trajectories are shown in a separate panel colored by the active regime.

## Usage

``` r
# S3 method for class 'pdmp_portrait'
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

  A `pdmp_portrait` object

- type:

  Character string selecting the plot variant (see Details).

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

The `type` argument selects the plot variant:

- `"default"`: a faceted layout with one panel per regime showing the
  ODE skeleton (vector field, nullclines, equilibria), plus a combined
  panel with all regime equilibria and switching trajectories.

- `"trajectories"`: only the switching trajectories, colored by the
  active regime, with switch points marked as dots.

- `"regime_NAME"`: single-regime portrait for the regime named NAME,
  showing its vector field, nullclines, and equilibria.

- `"combined"`: all regimes overlaid on a single panel. Each regime's
  nullclines get a distinct color. Use for visual comparison of regime
  dynamics.

Equilibrium points are encoded by both shape and color using the same
conventions as
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md):

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
telegraph <- system_spec(
    regimes = list(
        on  = list(x ~ k1 - d1 * x,  y ~ k2 * x - d2 * y),
        off = list(x ~ -d1 * x,       y ~ k2 * x - d2 * y)
    ),
    transitions = list(
        list(from = "on",  to = "off", rate = ~ lambda_off),
        list(from = "off", to = "on",  rate = ~ lambda_on)
    ),
    state_names = c("x", "y"),
    parms = list(k1 = 10, k2 = 5, d1 = 1, d2 = 0.5,
                 lambda_on = 0.5, lambda_off = 0.3),
    init = c(x = 5, y = 20),
    initial_regime = "on"
)

pp <- pdmp_portrait(telegraph, traj_length = 100)
plot(pp)
plot(pp, type = "trajectories")
plot(pp, type = "regime_on")
plot(pp, type = "combined")
} # }
```
