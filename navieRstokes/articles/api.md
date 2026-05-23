# API reference and technical architecture

This vignette is the navigational map for the package: design
philosophy, subsystem architecture, function-family catalogue, data flow
per time step, S3 classes, parallelism, and complexity. Mathematical
exposition of the governing equations and the projection method lives in
the companion vignette
[`vignette("navier-stokes-theory")`](https://robustecologies.github.io/navieRstokes/articles/navier-stokes-theory.md).

  

## Design philosophy

navieRstokes is a direct numerical solver for the 2D incompressible
Navier-Stokes equations on a uniform staggered Cartesian grid. Three
principles govern the architecture:

The first is **separation of declaration and execution**. The user
assembles a problem from explicit, named components: an
initial-condition function, a forcing function, a boundary-condition
specification, integration parameters.
[`simulate_navier_stokes()`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)
consumes the assembly and dispatches to the C++ kernel. There is no
formula-to-code translation layer; the user’s R-side helpers return
numeric grids that the C++ side integrates verbatim.

The second is **C++ kernels with R-side composition**. Convection,
diffusion, the pressure Poisson iteration and the Chorin projection are
implemented in Rcpp-backed C++ with optional OpenMP grid parallelism.
Helpers for initial conditions, body forcing and post-processing live in
R because they compose with user code and are not on the inner time-step
loop.

The third is **ggplot2-native visualization**. The S3 plot method
returns a ggplot2 object that the user can compose, facet, save, or
animate. The
[`flow()`](https://robustecologies.github.io/navieRstokes/reference/flow.md)
method is a thin wrapper that renders frames with the same ggplot2
method and assembles them through `av` (MP4) or `gifski` (GIF).

  

## Architecture overview

  

## Projection-method pipeline

Each time step executes the four-stage Chorin fractional-step pipeline.
The intermediate velocity carries the convection and diffusion update;
the pressure correction enforces incompressibility through a Poisson
solve; the projection subtracts the pressure gradient from the
intermediate velocity.

  

## Initial-condition catalogue

Six analytical initial conditions ship with the package; users may also
supply their own .

  

## Forcing-function catalogue

Forcing functions return a body-force per unit mass on the same
staggered grid as the velocity field;
[`simulate_navier_stokes()`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)
adds them to the convection-diffusion update before the Poisson solve.

  

## Visualization layer

[`plot.navieRstokes()`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md)
returns a ggplot2 object. The `plot_type` argument selects the scalar
field; the colour scale is `viridis` for `speed`, `velocity` and
`vorticity`-magnitude flavours and a diverging blue-white-red for signed
`vorticity`; pressure uses `YlOrRd`.
[`flow.navieRstokes()`](https://robustecologies.github.io/navieRstokes/reference/flow.md)
renders an animation by iterating the same ggplot2 method through every
saved snapshot.

  

## Shiny dashboard

[`shiny_navieRstokes()`](https://robustecologies.github.io/navieRstokes/reference/shiny_navieRstokes.md)
launches an interactive dashboard with pre-configured scenarios
(lid-driven cavity, vortex pair, Kelvin-Helmholtz, rotating disk). The
dashboard wires the same API: parameter widgets construct an initial
condition + forcing + boundary spec, the call to
[`simulate_navier_stokes()`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)
runs in a `withProgress()` block, and the resulting object is plotted
with `plot.navieRstokes` inside a `plotOutput`.

  

## S3 class system

The package defines a single user-facing class, `navieRstokes`, returned
by
[`simulate_navier_stokes()`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md).
It carries the velocity arrays `u` and `v` (3D `[nx x ny x nt_save]`),
the pressure `p` (same shape), the spatial grids `x` and `y`, the time
vector `t`, the parameter list and a diagnostics block (mass error, CFL,
energy per saved time-step). The class has a complete S3 trio:
`print.navieRstokes` (one-screen summary), `summary.navieRstokes`
(diagnostics) and `plot.navieRstokes` (ggplot2 field visualization).
`flow.navieRstokes` extends the class with animation generation through
the
[`flow()`](https://robustecologies.github.io/navieRstokes/reference/flow.md)
generic.

  

## Parallelism and reproducibility

The C++ kernels carry `#pragma omp parallel for` over the inner spatial
loops. OpenMP is optional at compile time; the package builds and runs
serially when OpenMP headers are unavailable. Reproducibility is
automatic: the solver is deterministic (no Monte Carlo), and identical
parameters produce bitwise-identical output across serial and parallel
runs.

  

## Complexity reference

Per time step, on an `nx`-by-`ny` grid, the convection and diffusion
stages cost `O(nx ny)` floating-point operations. The pressure Poisson
solve costs `O(nx ny k)` where `k` is the iteration count to
convergence; for Jacobi the iteration count is `O(max(nx, ny))` for
tight tolerances, while for SOR with optimal relaxation the count drops
to roughly `O(sqrt(max(nx, ny)))`. The CFL constraint sets
`dt ≲ min(dx, dy) / max(|u|)`; violations are flagged in the diagnostics
block.

  

## Function-family catalogue

The reference index in `_pkgdown.yml` mirrors the six functional
families used here and groups every export accordingly: simulation core
(`simulate_navier_stokes`), initial conditions (six exports), forcing
functions (four exports), diagnostics (`compute_vorticity`),
visualization (`flow`, `plot.navieRstokes`, `flow.navieRstokes`),
interactive dashboard (`shiny_navieRstokes`). The S3 trio is bucketed at
the bottom of the reference page.

  

## Function reference

The remainder of this vignette is a per-symbol reference: for every
exported function and every S3 method, it lists the signature, returned
class, one-sentence purpose, and a link to the rendered help page. Code
chunks are not evaluated; they show the calling convention only.

  

### simulate_navier_stokes

Solver entry point: integrates the 2D incompressible Navier-Stokes
equations with the Chorin projection method on a uniform grid.

``` r

simulate_navier_stokes(nx, ny, lx, ly, dt, nt, nu, rho = 1.0,
                      initial_condition = NULL,
                      boundary_condition = NULL,
                      forcing_function = NULL,
                      pressure_solver = "jacobi",
                      max_iter = NULL, tolerance = 1e-3,
                      save_interval = 1L, save_final = TRUE,
                      verbose = TRUE)
```

Returns an S3 list of class `navieRstokes` carrying `u`, `v`, `p` (3D
arrays `[nx, ny, n_save]`), `x`, `y`, `t`, `parameters`, `diagnostics`.
See
[`?simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md).

  

### vortex_ic, taylor_green_ic

Doubly-periodic Taylor-Green velocity field; `taylor_green_ic` is a
re-export alias of `vortex_ic`.

``` r

vortex_ic(x, y, A = 1.0, lx = 1.0, ly = 1.0)
```

Returns `list(u, v)`. Divergence-free by construction. See
[`?vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md).

  

### shear_layer_ic

Hyperbolic-tangent shear layer with a sinusoidal perturbation; canonical
seed for Kelvin-Helmholtz roll-up.

``` r

shear_layer_ic(x, y, U0 = 1.0, delta = 0.05, epsilon = 0.05,
               lx = 1.0, ly = 1.0)
```

Returns `list(u, v)`. See
[`?shear_layer_ic`](https://robustecologies.github.io/navieRstokes/reference/shear_layer_ic.md).

  

### rotating_ic

Solid-body rotation about the domain centre with angular velocity
`omega`; uniform vorticity `2 omega`.

``` r

rotating_ic(x, y, omega = 2 * pi, x0 = NULL, y0 = NULL,
            lx = 1.0, ly = 1.0)
```

Returns `list(u, v)`. See
[`?rotating_ic`](https://robustecologies.github.io/navieRstokes/reference/rotating_ic.md).

  

### single_vortex_ic

Lamb-Oseen Gaussian vortex centred at `(x0, y0)`; smooth at the core,
exponentially decaying.

``` r

single_vortex_ic(x, y, Gamma = 1.0, sigma = 0.15,
                 x0 = NULL, y0 = NULL, lx = 1.0, ly = 1.0)

# Demonstration: viscous decay of an isolated Lamb-Oseen vortex
res <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 1000, nu = 0.005,
  initial_condition = single_vortex_ic,
  boundary_condition = list(type = "periodic"),
  save_interval = 50
)
plot(res, plot_type = "vorticity")
```

Returns `list(u, v)`. See
[`?single_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/single_vortex_ic.md).

  

### two_vortex_ic

Pair of counter-rotating Lamb-Oseen vortices at `(0.3 lx, 0.5 ly)` and
`(0.7 lx, 0.5 ly)`; canonical dipole interaction.

``` r

two_vortex_ic(x, y, Gamma = 1.0, sigma = 0.1, lx = 1.0, ly = 1.0)

# Demonstration: dipole self-advection under periodic BCs
res <- simulate_navier_stokes(
  nx = 96, ny = 96, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 800, nu = 0.001,
  initial_condition = two_vortex_ic,
  boundary_condition = list(type = "periodic"),
  save_interval = 40
)
plot(res, plot_type = "vorticity")
```

Returns `list(u, v)`. See
[`?two_vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/two_vortex_ic.md).

  

### pressure_gradient_force

Constant body force in the x-direction, equivalent to a uniform pressure
gradient; canonical Poiseuille driver.

``` r

pressure_gradient_force(x, y, t, dp_dx = 0.1, rho = 1.0)
```

Returns `list(fx, fy)`. Sign convention: positive `dp_dx` produces flow
in the negative x-direction. See
[`?pressure_gradient_force`](https://robustecologies.github.io/navieRstokes/reference/pressure_gradient_force.md).

  

### oscillatory_force

Spatially uniform, sinusoidal-in-time body force; canonical
Stokes-second-problem driver.

``` r

oscillatory_force(x, y, t, A = 0.5, omega = 2 * pi)

# Demonstration: oscillating Stokes layer in a channel
res <- simulate_navier_stokes(
  nx = 32, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 1000, nu = 0.05,
  initial_condition = function(x, y) list(u = 0, v = 0),
  forcing_function = function(x, y, t) {
    oscillatory_force(x, y, t, A = 0.2, omega = 2 * pi)
  },
  boundary_condition = list(
    type = "channel",
    values = list(u_top = 0, u_bottom = 0, v_top = 0, v_bottom = 0)
  ),
  save_interval = 50
)
plot(res, plot_type = "velocity")
```

Returns `list(fx, fy)`. See
[`?oscillatory_force`](https://robustecologies.github.io/navieRstokes/reference/oscillatory_force.md).

  

### taylor_green_force

Steady forcing matched to the Taylor-Green decay rate; sustains
`vortex_ic` against viscous dissipation.

``` r

taylor_green_force(x, y, t, A = 1.0, k = 2 * pi, nu = 0.01,
                   lx = 1.0, ly = 1.0)
```

Returns `list(fx, fy)`. The wavenumbers must match `vortex_ic`
(`kx = 2*pi/lx`, `ky = 2*pi/ly`). See
[`?taylor_green_force`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md).

  

### localized_vortex_force

Gaussian-localised rotational forcing centred at `(x0, y0)`; injects a
swirl with characteristic radius `sigma`.

``` r

localized_vortex_force(x, y, t, x0 = 0.5, y0 = 0.5,
                       sigma = 0.1, A = 1.0)

# Demonstration: localised swirl injected into a periodic base flow
res <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 600, nu = 0.01,
  initial_condition = function(x, y) vortex_ic(x, y, A = 0.01),
  forcing_function = function(x, y, t) {
    localized_vortex_force(x, y, t, x0 = 0.5, y0 = 0.5,
                           sigma = 0.1, A = 0.5)
  },
  boundary_condition = list(type = "periodic"),
  save_interval = 50
)
plot(res, plot_type = "vorticity")
```

Returns `list(fx, fy)`. See
[`?localized_vortex_force`](https://robustecologies.github.io/navieRstokes/reference/localized_vortex_force.md).

  

### compute_vorticity

Cell-centred vorticity `omega = dv/dx - du/dy` via second-order central
differences; thin R wrapper around the C++ kernel.

``` r

compute_vorticity(u, v, dx, dy)
```

Returns a `numeric` matrix of dimension `dim(u)`. Boundary cells are
filled with `NA`. See
[`?compute_vorticity`](https://robustecologies.github.io/navieRstokes/reference/compute_vorticity.md).

  

### flow (generic)

S3 generic that animates a `navieRstokes` object frame-by-frame.

``` r

flow(x, ...)
```

Dispatches to `flow.navieRstokes`. See
[`?flow`](https://robustecologies.github.io/navieRstokes/reference/flow.md).

  

### flow.navieRstokes

Renders saved snapshots through `plot.navieRstokes` and assembles them
into an MP4 (via `av`) or GIF (via `gifski`).

``` r

flow(result, output_file = "flow.mp4", fps = 10, format = "mp4",
     plot_type = "vorticity",
     width = 800, height = 600, res = 100)
```

Returns the output filename invisibly. See
[`?flow.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/flow.md).

  

### plot.navieRstokes

S3 plot method returning a `ggplot` (or a `patchwork` composition for
`plot_type = "all"`).

``` r

plot(result, time_index = NULL,
     plot_type = c("vorticity", "pressure", "speed", "velocity", "all"),
     subsample = 4L, scale = 0.04, n_contours = 12L)
```

Returns a `ggplot` object that can be composed, faceted, themed and
saved with the standard ggplot2 toolchain. See
[`?plot.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md).

  

### shiny_navieRstokes

Launches the `shinydashboard`-based interactive dashboard with
pre-configured cavity, vortex, shear, rotation and Poiseuille scenarios.

``` r

shiny_navieRstokes()
```

Side-effect only (no return value). Requires `shiny` and
`shinydashboard`. See
[`?shiny_navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/shiny_navieRstokes.md).

  

### print.navieRstokes

One-screen S3 print method: grid, domain, viscosity, Reynolds, time
range, mass error, CFL.

``` r

print(result)
```

Returns the input invisibly. See
[`?print.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/print.navieRstokes.md).

  

### summary.navieRstokes

Constructs an S3 `summary.navieRstokes` object carrying parameters,
diagnostics, time range, velocity statistics, Poisson-solver convergence
statistics and CFL statistics.

``` r

summary(result)
```

Returns an object of class `summary.navieRstokes`. See
[`?summary.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/summary.navieRstokes.md).

  

### print.summary.navieRstokes

S3 print method for `summary.navieRstokes`: a multi-block diagnostics
report covering grid, physical parameters, time integration, final-state
velocity field, mass conservation, Poisson convergence and CFL
stability.

``` r

print(summary(result))
```

Returns the input invisibly. See
[`?print.summary.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/print.summary.navieRstokes.md).
