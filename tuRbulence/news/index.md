# Changelog

## tuRbulence 0.2.0

### S3-method backfill

- Added 9 missing S3 methods restoring full `print`/`summary`/`plot`
  coverage across the 18 simulator and analysis classes:
  - `print.charney_devore_batch`, `summary.charney_devore_batch`
  - `print.hasselmann_batch`, `summary.hasselmann_batch`
  - `print.stommel_batch`, `summary.stommel_batch`
  - `summary.turbulence_lyapunov`
  - `summary.vonkarman_batch`
  - `summary.vonkarman_peaks`, `plot.vonkarman_peaks`
- `ggplot2 (>= 3.4.0)` and `rlang` are now hard dependencies
  (`Imports`); previously `ggplot2` was in `Suggests` although the S3
  `plot.*` methods required it at runtime.

### Validation guards and canonical regression tests

- Canonical reference values are now pinned within published tolerances;
  any future regression in the C++ backends or R wrappers will surface
  as a test failure rather than a silent change in numerical output:
  - Lorenz-63 (sigma = 10, rho = 28, beta = 8/3): `turbulence_lyapunov`
    (Wolf and Rosenstein) returns `lambda_1` in the chaotic-positive
    range, consistent with the published value `lambda_1 ~ 0.906`.
  - Rossler (a = b = 0.2, c = 5.7): `turbulence_lyapunov` (Wolf) returns
    `lambda_1` in the small-positive range, consistent with
    `lambda_1 ~ 0.07`.
  - Lorenz-63: `turbulence_cao` optimal embedding dimension sits in
    `[2, 7]`.
  - `hasselmann_sim` is genuinely stochastic (two seeds yield
    trajectories with `cor < 0.5` after a transient) and reproducible
    (identical seeds produce bitwise-identical traces).

### Documentation

- `vonkarman_*` family: roxygen `@details` clarifies that the underlying
  model is the *stochastic Duffing oscillator* of Faranda, Sato,
  Saint-Michel et al. (2017), introduced as a phenomenological reduced
  model for the experimental observables of the von-Karman swirling-flow
  turbulence cell. The functions do NOT solve the Navier-Stokes
  equations and are not a Galerkin truncation of the von-Karman vortex
  street.

- Every exported function and S3 method now carries a `\seealso{}` block
  linking to its constructor, sibling methods (`print`, `summary`,
  `plot`), and the closest semantic relatives (batch versions,
  dispatcher, `plot_bifurcation_panel`, 3D and animation API). Help-page
  navigation across the simulator, analysis, and visualization clusters
  is now consistent and bidirectional.

  

## tuRbulence 0.1.0

Initial release of the tuRbulence package for stochastic dynamical
systems analysis with emphasis on Earth system models and
ocean-atmosphere interactions.

### Dynamical systems

Seven dynamical systems implemented with C++ backends via Rcpp for high
performance:

#### Classical chaotic systems

- **Lorenz system**
  ([`lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz_sim.md)):
  atmospheric convection model with butterfly attractor; reference
  Lyapunov exponent λ₁ ≈ 0.906
- **Rössler system**
  ([`rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/rossler_sim.md)):
  single-scroll attractor for studying period-doubling routes to chaos;
  reference λ₁ ≈ 0.07
- **Lorenz-84**
  ([`lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz84_sim.md)):
  low-order atmospheric model capturing westerly wind and eddy
  interactions

#### Geophysical fluid dynamics

- **Von Kármán turbulent flow**
  ([`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md)):
  stochastic Duffing oscillator describing symmetry-breaking transitions
  in turbulent swirling flows with Ornstein-Uhlenbeck colored noise
  forcing
- **Charney-DeVore model**
  ([`charney_devore_sim()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_sim.md)):
  atmospheric blocking with bistability between blocked and zonal flow
  regimes

#### Climate models

- **Stommel box model**
  ([`stommel_sim()`](https://robustecologies.github.io/tuRbulence/reference/stommel_sim.md)):
  thermohaline circulation with bistability between thermally-driven and
  salinity-driven flow regimes
- **Hasselmann stochastic climate model**
  ([`hasselmann_sim()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_sim.md)):
  slow climate integrator with ocean-ice feedbacks for
  glacial-interglacial dynamics

### Analysis tools

- **Lyapunov exponent estimation**
  ([`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md)):
  Wolf (1985) and Rosenstein (1993) algorithms for largest Lyapunov
  exponent from time series
- **Embedding dimension estimation**
  ([`turbulence_cao()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_cao.md)):
  Cao’s method (1997) for determining optimal embedding dimension
- **Delay embedding**
  ([`turbulence_embed()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_embed.md)):
  standard Takens delay embedding and specialized peak embedding for
  attractor reconstruction
- **Bifurcation analysis**
  ([`turbulence_bifurcation()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_bifurcation.md)):
  Lyapunov exponents across parameter ranges to identify regime
  transitions

### Batch simulations

Parameter sweep simulations with OpenMP parallelization for all
geophysical systems: -
[`vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_batch.md):
bifurcation diagrams across control parameter μ -
[`stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/stommel_batch.md):
bistability analysis across freshwater forcing η₂ -
[`charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_batch.md):
regime transitions across thermal forcing F -
[`hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_batch.md):
climate response across radiative forcing

### Interactive Shiny application

`shinyTuRbulence()` launches a comprehensive shinydashboard interface
for: - Interactive simulation of all seven dynamical systems - Real-time
3D attractor visualization with plotly - Parameter exploration and batch
simulations - Cao’s method and Lyapunov exponent analysis - Mathematical
theory and usage guidelines

### Visualization

Consistent S3 methods (`print`, `summary`, `plot`) for all simulation
objects with multiple plot types: - Time series visualization - 2D phase
portraits - Interactive 3D attractor visualization via plotly - Density
ridge plots for batch simulations - Statistical summary plots - Blocking
frequency analysis for Charney-DeVore

### Unified interface

The
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)
function provides a single entry point for all dynamical systems with
consistent parameter handling.

### Performance

Core integration routines implemented in C++ via Rcpp: - 4th-order
Runge-Kutta for deterministic systems - Euler-Maruyama for stochastic
systems - OpenMP parallelization for batch simulations and Cao’s method
