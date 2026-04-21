# Visualise a simulated dynamical system

Produces ggplot2 or plotly visualisations of a simulation produced by
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md).
The `type` argument selects the view (time series, 2D phase portrait, 3D
attractor, or PDE heatmap / snapshot) and the method dispatches to the
appropriate rendering backend. An optional sensitivity overlay
re-simulates the system with a perturbed initial condition and draws the
perturbed trajectory in a contrasting colour, providing a visual
diagnostic of sensitive dependence on initial conditions.

## Usage

``` r
# S3 method for class 'dyn_sim'
plot(x, type = NULL, show_transient = NULL, times = NULL, title = NULL, ...)
```

## Arguments

- x:

  A `dyn_sim` object produced by
  [`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md).

- type:

  Character; plot type. Admissible values depend on the model class (see
  Details). Defaults to `"timeseries"` for non-spatial models and
  `"heatmap"` for PDE/RDME models.

- show_transient:

  Logical; if `TRUE` the full trajectory including transients is used,
  otherwise only the post-transient attractor. Defaults to `TRUE` for
  `type = "timeseries"` and `FALSE` otherwise.

- times:

  Numeric vector of times at which to render snapshots for
  `type = "snapshot"` (PDE only). If `NULL`, four evenly spaced values
  covering the simulated interval are selected.

- title:

  Character plot title. If `NULL` (default), the model name stored in
  `x$model$meta$name` is used; otherwise a generic fallback title is
  provided.

- ...:

  Additional arguments forwarded to the backend plot function, most
  notably `attractor_mode` which selects `"lines"` or `"markers"`
  rendering for the 3D attractor (default: `"markers"` for stochastic
  solvers, `"lines"` otherwise).

- epsilon:

  Numeric scalar, named numeric vector, or `NULL` (default). Controls
  the perturbed-initial-condition overlay; see Details.

- vars:

  Character vector of length 2 naming the two state variables projected
  by `type = "phase"`. Defaults to the first two state variables of the
  model. An error is raised if either name is not a state of `x$model`.

- time_colored:

  Logical used by `type = "attractor"`; if `TRUE` (default) the 3D
  trajectory is coloured by time with Viridis and a colourbar is shown.
  When `epsilon` is supplied this is forced to `FALSE`, so that the
  original and perturbed trajectories are distinguished by colour.

- perturbed_color:

  Character colour used for the perturbed trajectory when `epsilon` is
  non-`NULL`. Default `"#d9c22b"` (dark goldenrod). The original
  trajectory is drawn in `"#2a1766"` (deep purple).

## Value

For ggplot2 backends (time series, phase, PDE heatmap, PDE snapshot,
PDMP time series, spatially averaged PDE time series) a `ggplot` object
is returned invisibly after being printed. For the plotly backend (3D
attractor) a `plotly` object is returned directly (and will self-render
at the R console).

## Details

Visualise a simulated dynamical system

The set of admissible plot types depends on the nature of the model.

For ODE, SDE, DDE, map, SSA/tau-leap, jump-diffusion and PDMP
simulations:

- `"timeseries"`:

  (default) State variables against time (or iteration number, for
  maps), one colour per variable. PDMP trajectories are coloured by
  active regime instead.

- `"phase"`:

  Two-dimensional phase-plane projection on a pair of state variables
  selected by `vars`. Rendered with ggplot2.

- `"attractor"`:

  Three-dimensional attractor projection on the first three state
  variables, rendered with plotly. When `time_colored = TRUE` the
  trajectory is colour-coded by time with the Viridis colour map.

For 1D PDE, 2D PDE and RDME models:

- `"heatmap"`:

  (default) Space-time heatmap for 1D problems; spatial heatmap at the
  final time for 2D problems and lattice RDME.

- `"snapshot"`:

  Line plots of the spatial profile at the times given by `times` (four
  evenly spaced values by default).

- `"timeseries"`:

  Spatially averaged state variables plotted against time.

**Sensitivity overlay.** The overlay is triggered by setting `epsilon`
to a non-`NULL` value. The package re-invokes
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
with the initial state perturbed by `epsilon`, preserving solver,
parameters, `t_max`, `discard_transient` and random seed, then draws the
perturbed trajectory on top of the original with a contrasting colour
and an explanatory subtitle and legend. `epsilon` accepts two forms:

- a numeric scalar \\\varepsilon\\: every state is shifted by the same
  amount, \\y_0 \leftarrow y_0 + \varepsilon \\ \mathbf{1}\\;

- a named numeric vector indexed by state names: only the named
  components are shifted, e.g. `c(X = 0.001, Z = -0.01)`.

The overlay is not supported for PDE, RDME and PDMP plots; passing
`epsilon` in those cases emits a warning and the argument is ignored.

## References

Strogatz, S. H. (2015). *Nonlinear Dynamics and Chaos*, 2nd edition.
Westview Press. ISBN 9780813349107.

Hastings, A. and Powell, T. (1991). Chaos in a three-species food chain.
*Ecology*, 72(3), 896-903.
[doi:10.2307/1940591](https://doi.org/10.2307/1940591) .

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*, 2nd
edition. Springer.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)
.

## See also

[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
for the constructor that produces the input object;
[`print.dyn_sim`](https://robustecologies.github.io/janos/reference/print.dyn_sim.md)
and
[`summary.dyn_sim`](https://robustecologies.github.io/janos/reference/summary.dyn_sim.md)
for the companion S3 methods;
[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md)
for defining the underlying dynamical system;
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md),
[`map_portrait`](https://robustecologies.github.io/janos/reference/map_portrait.md),
[`sde_portrait`](https://robustecologies.github.io/janos/reference/sde_portrait.md),
[`dde_portrait`](https://robustecologies.github.io/janos/reference/dde_portrait.md)
and
[`pdmp_portrait`](https://robustecologies.github.io/janos/reference/pdmp_portrait.md)
for vector-field based portraits that complement the trajectory-based
views provided here.

## Examples

``` r
if (FALSE) { # \dontrun{
# Hastings-Powell tri-trophic food chain (chaotic "teacup" attractor)
hp <- model_spec(
    rhs = list(
        X ~ X * (1 - X) - (a1 * X * Y) / (1 + b1 * X),
        Y ~ (a1 * X * Y) / (1 + b1 * X) -
            (a2 * Y * Z) / (1 + b2 * Y) - d1 * Y,
        Z ~ (a2 * Y * Z) / (1 + b2 * Y) - d2 * Z
    ),
    state_names = c("X", "Y", "Z"),
    parms = list(a1 = 5, b1 = 3, a2 = 0.1, b2 = 2,
                 d1 = 0.4, d2 = 0.01),
    init  = c(X = 0.7, Y = 0.2, Z = 9.0)
)
run <- dyn_sim(hp, t_max = 1000, solver = solver_rk45(),
               discard_transient = 600, verbose = FALSE)

# time series of every state variable
plot(run, type = "timeseries")

# 2D phase-plane projection on the (X, Z) pair
plot(run, type = "phase", vars = c("X", "Z"))

# 3D attractor, colour-coded by time
plot(run, type = "attractor")

# 3D attractor with uniform colouring
plot(run, type = "attractor", time_colored = FALSE)

# sensitivity overlay: shift every state by 10^-3
plot(run, type = "attractor", epsilon = 1e-3)

# sensitivity overlay: perturb only X
plot(run, type = "phase", vars = c("X", "Z"),
     epsilon = c(X = 1e-3))
} # }
```
