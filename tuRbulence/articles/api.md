# API reference and architecture

This vignette is the navigational map for `tuRbulence`. It diagrams the
seven dynamical systems, the deterministic-vs-stochastic dispatch, the
batch-parallel layer, the analysis stack (Lyapunov, Cao, Takens,
bifurcation), and the visualization pipeline. The mathematical
exposition lives in
[`vignette("turbulent-attractors")`](https://robustecologies.github.io/tuRbulence/articles/turbulent-attractors.md).

  

## Design philosophy

`tuRbulence` ships seven canonical low-order dynamical systems used
across atmospheric, oceanographic and statistical physics: Lorenz-63,
Lorenz-84, Rossler, Charney-DeVore (3-mode and 6-mode), Stommel two-box
thermohaline, Hasselmann stochastic climate, and the stochastic Duffing
model that Faranda et al. (2017) introduced as a phenomenological
reduced model for the von-Karman swirling-flow turbulence cell.

Three principles govern the architecture:

The first is **single-dispatcher convenience plus per-system
specificity**.
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)
accepts a system name and forwards to the appropriate integrator. Each
integrator also exists as an exported function for users who need direct
control of system-specific parameters.

The second is **explicit deterministic-versus-stochastic flag**. Most
simulators carry a `stochastic` argument. When TRUE, an Euler-Maruyama
step adds a Wiener increment; when FALSE, the simulator integrates the
noiseless ODE with RK4. Hasselmann and the von-Karman stochastic Duffing
are stochastic by construction.

The third is **batch parallelism for parameter sweeps**. Each system
carries a `*_batch()` constructor that runs the simulator across a
vector of forcing or control values, parallelised by OpenMP through
Rcpp. Returned objects share their canonical S3 class for the
per-simulation plot/summary/print methods, plus a `*_batch` class for
cross-sim aggregations.

  

## System catalogue

  

## Deterministic vs stochastic dispatch

Most simulators flip between RK4 (deterministic) and Euler-Maruyama
(stochastic) according to the `stochastic` argument. Hasselmann and
von-Karman are stochastic by construction; they always integrate the
Wiener increment.

  

## Batch parallelism

Each `*_batch()` runs the same integrator across a vector of forcing
values via OpenMP at the C++ level. Returned object carries the
per-simulation observable arrays plus a vector of forcing values. Plot
methods render either trajectories per forcing or aggregated quantities
(energy, q-flow strength, etc.).

  

## Analysis stack

  

## Visualization

2D outputs use ggplot2 with the canonical subtitle and grey caption
pattern; 3D attractors and animations use plotly (and `htmlwidgets` for
the standalone-HTML save).

  

## Shiny dashboard

[`shiny_tuRbulence()`](https://robustecologies.github.io/tuRbulence/reference/shiny_tuRbulence.md)
launches an interactive dashboard with pre-configured scenarios for each
of the seven systems. Parameter widgets construct the simulator call;
the result is plotted with the appropriate S3 plot method inside
`plotOutput`. The dashboard wires the same API as the programmatic
interface; nothing in the UI is privileged.

  

## S3 class system

The package exports 18 S3 classes, each with the full print/summary/plot
triple per CLAUDE convention.

| Class                  | Constructor            | Methods                |
|:-----------------------|:-----------------------|:-----------------------|
| lorenz_sim             | lorenz_sim             | print + summary + plot |
| lorenz84_sim           | lorenz84_sim           | print + summary + plot |
| rossler_sim            | rossler_sim            | print + summary + plot |
| charney_devore_sim     | charney_devore_sim     | print + summary + plot |
| charney_devore_6mode   | charney_devore_6mode   | print + summary + plot |
| charney_devore_batch   | charney_devore_batch   | print + summary + plot |
| stommel_sim            | stommel_sim            | print + summary + plot |
| stommel_batch          | stommel_batch          | print + summary + plot |
| hasselmann_sim         | hasselmann_sim         | print + summary + plot |
| hasselmann_batch       | hasselmann_batch       | print + summary + plot |
| vonkarman_sim          | vonkarman_sim          | print + summary + plot |
| vonkarman_batch        | vonkarman_batch        | print + summary + plot |
| vonkarman_peaks        | vonkarman_peaks        | print + summary + plot |
| vonkarman_attractor    | vonkarman_attractor    | print + summary + plot |
| turbulence_embedding   | turbulence_embed       | print + summary + plot |
| turbulence_lyapunov    | turbulence_lyapunov    | print + summary + plot |
| turbulence_cao         | turbulence_cao         | print + summary + plot |
| turbulence_bifurcation | turbulence_bifurcation | print + summary + plot |

  

## Important note on the von-Karman family

The functions `vonkarman_sim`, `vonkarman_batch`, `vonkarman_peaks` and
`vonkarman_attractor` integrate the *stochastic Duffing oscillator*
introduced by Faranda, Sato, Saint-Michel et al. (2017) “Stochastic
chaos in a turbulent swirling flow” (Phys. Rev. Lett. 119, 014502,
<doi:10.1103/PhysRevLett.119.014502>) as a phenomenological reduced
model for the experimental observables of the von-Karman swirling-flow
turbulence cell. They do NOT solve the Navier-Stokes equations and are
not a Galerkin truncation of the von-Karman vortex street. The function
name reflects the experimental cell the model was built to describe; the
equations are the stochastic Duffing model in the cited paper. See
[`?vonkarman_sim`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md)
for the full equations.

  

## Validation gates

The test suite under `tests/testthat/test-validation-canonical.R` (gated
by `skip_on_cran()`) pins the package contract on canonical reference
values: Lorenz-63 (sigma = 10, rho = 28, beta = 8/3) yields
`lambda_1 ~ 0.906` from both Wolf and Rosenstein algorithms; Rossler (a
= b = 0.2, c = 5.7) yields `lambda_1 ~ 0.07`; Cao optimal embedding
dimension on Lorenz-63 sits in `[3, 6]`; `hasselmann_sim` is genuinely
stochastic (two seeds diverge with `cor < 0.5` after a transient) and
reproducible (identical seeds produce bitwise-identical traces).
