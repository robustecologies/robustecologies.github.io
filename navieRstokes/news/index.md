# Changelog

## navieRstokes 0.3.2

### CRAN-readiness fixes

- Package title shortened to comply with the CRAN 65-character limit;
  sentence-case form retained.

- `LazyData` field set to `false` (the package ships no datasets),
  removing a CRAN NOTE.

  

### Documentation

- Help pages now expose a complete navigation cluster via `@seealso`:
  every initial-condition function cross-references its siblings and the
  solver; every forcing function cross-references its siblings and the
  solver; the S3 trio (`print.navieRstokes`, `summary.navieRstokes`,
  `plot.navieRstokes`, `flow.navieRstokes`) cross-references each other
  and the constructor; `compute_vorticity` joins the visualisation
  cluster.

  

## navieRstokes 0.3.1

### Visualization improvements

- `plot.navieRstokes` gains `plot_type = "all"` for a three-panel
  patchwork composition showing vorticity, speed and pressure
  simultaneously. `patchwork` is now a hard dependency (`Imports`).

- Contour lines (`geom_contour`) are overlaid on vorticity and pressure
  plots by default. New argument `n_contours` (default 12) controls the
  number of iso-contour levels. Set to `0L` to suppress.

- Raster cells are now rendered uninterpolated (`interpolate = FALSE`),
  giving each grid cell its exact discrete colour value and restoring
  the granularity of the previous
  [`image()`](https://rdrr.io/r/graphics/image.html)-based output.

  

### Bug fixes

- Centerline velocity profiles in `vignettes/lid-driven-cavity.Rmd` were
  extracting a horizontal slice instead of the expected vertical
  centerline, causing zig-zag line artifacts. The extraction is now
  `u_final[i_mid, ]` (all y at fixed x) and `v_final[, j_mid]` (all x at
  fixed y), consistent with the Ghia et al. benchmark convention.

- Poiseuille flow analytical solution in
  `vignettes/advanced-examples.Rmd` referenced `parameters$rho`, which
  is not stored by the incompressible solver. The formula is now
  `G / (2 * nu) * y * (h - y)` using kinematic viscosity only, matching
  the unit-density incompressible NS formulation.

  

## navieRstokes 0.3.0

### Visualization migrated to ggplot2

- **BREAKING:** `plot.navieRstokes` and `flow.navieRstokes` now build
  their output through `ggplot2` instead of base graphics. The plot
  method returns a `ggplot` object that the user can compose, facet or
  save with any ggplot2-aware tool. The `plot_type` argument is widened
  to `"vorticity"`, `"pressure"`, `"speed"`, `"velocity"`; the velocity
  flavour overlays a quiver of `geom_segment` arrows on the speed
  raster. New `subsample` and `scale` arguments control quiver density
  and arrow length.

- `flow.navieRstokes` switches to the same ggplot2 frame renderer and
  prefers
  [`ragg::agg_png()`](https://ragg.r-lib.org/reference/agg_png.html)
  when the `ragg` package is installed, eliminating `mbcsToSbcs`
  warnings on Greek labels under common locales.

- `ggplot2 (>= 3.4.0)` and `rlang` are now hard dependencies
  (`Imports`); `grDevices` and `graphics` are no longer imported. Saved
  figures requiring Greek labels should use
  [`ragg::agg_png`](https://ragg.r-lib.org/reference/agg_png.html)/[`ragg::agg_tiff`](https://ragg.r-lib.org/reference/agg_tiff.html).

  

## navieRstokes 0.2.0

### Major new feature: OpenMP parallel computation

This release introduces OpenMP parallelization for all major
computational kernels, enabling significant performance improvements on
multi-core systems for large grids.

#### OpenMP implementation details

- **Parallelized kernels**: All computationally intensive C++ kernels
  now support OpenMP parallelization:

  - `solve_poisson_jacobi_cpp()`: Jacobi iterations for pressure Poisson
    equation
  - `compute_convection_cpp()`: Upwind scheme for convective terms
  - `compute_diffusion_cpp()`: Laplacian operator for diffusion
  - `compute_divergence_cpp()`: Velocity field divergence
  - `compute_gradient_cpp()`: Pressure gradient computation
  - `compute_vorticity_cpp()`: Vorticity calculation

- **Conditional parallelization**: OpenMP overhead can exceed benefits
  for small grids. The package implements intelligent threshold-based
  activation: parallelization is only enabled for grids with \>= 16,384
  cells (128×128 or larger). Smaller grids execute serially for optimal
  performance.

- **Loop collapsing**: All 2D nested loops use `collapse(2)` directive
  for better workload distribution across threads.

- **Memory optimization**: Jacobi solver uses pre-allocated
  double-buffering instead of per-iteration memory allocation, reducing
  overhead significantly.

#### Performance characteristics

| Grid size | Cells | Serial time (1000 steps) | Notes |
|----|----|----|----|
| 64×64 | 4,096 | ~0.4-0.7s | Below threshold, serial execution |
| 128×128 | 16,384 | ~2.5s | OpenMP threshold, parallelization active |
| 256×256 | 65,536 | ~30s | Full OpenMP benefit |

**Important**: OpenMP benefit depends on problem characteristics. For
well-conditioned problems where the Poisson solver converges quickly
(1-4 iterations), the parallelization overhead may not be fully
amortized. Maximum benefit occurs for: - Large grids (256×256 or
larger) - Problems requiring many Poisson iterations - Systems with
multiple physical cores

#### System requirements

OpenMP is optional and automatically detected at compile time. The
package compiles and runs correctly without OpenMP support, falling back
to serial execution.

- **Linux/macOS**: Requires OpenMP-capable compiler (GCC 4.2+, Clang
  with libomp)
- **Windows**: Rtools includes OpenMP support by default

### Bug fixes and improvements

#### Critical fixes

- **Taylor-Green IC/forcing consistency**: Fixed sign convention
  mismatch between
  [`vortex_ic()`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)
  and
  [`taylor_green_force()`](https://robustecologies.github.io/navieRstokes/reference/taylor_green_force.md).
  Both now use the same convention: u = A*sin(kx)*cos(ky), v =
  -A*cos(kx)*sin(ky).

- **Poisson solver compatibility conditions**: Added mean-zero RHS
  enforcement for Neumann and periodic boundary conditions (mathematical
  compatibility condition for Poisson equation with pure Neumann BCs).

- **Boundary conditions at t=0**: Fixed boundary conditions not being
  applied to initial state, which caused incorrect initial conditions
  for lid-driven cavity simulations.

#### Enhancements

- **Domain scaling for initial conditions**: Added `lx` and `ly`
  parameters to
  [`vortex_ic()`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
  [`rotating_ic()`](https://robustecologies.github.io/navieRstokes/reference/rotating_ic.md),
  and
  [`shear_layer_ic()`](https://robustecologies.github.io/navieRstokes/reference/shear_layer_ic.md)
  for non-unit domains.

- **[`taylor_green_ic()`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)
  alias**: Added alias for
  [`vortex_ic()`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md)
  for clearer naming.

- **Residual-based Poisson convergence**: Changed from update-based to
  residual-based convergence criterion for more physically meaningful
  convergence monitoring.

- **Periodic BC corner handling**: Explicit corner handling for periodic
  boundary conditions to avoid order-dependent overwrites.

- **NA at boundaries for divergence**: `compute_divergence_cpp()` now
  sets boundary values to NA instead of zero for more accurate RMS
  diagnostics.

- **CFL guard for nu=0**: Added protection against division by zero in
  diffusion CFL calculation when viscosity is zero.

- **`save_final` parameter**: New parameter in
  [`simulate_navier_stokes()`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)
  to control whether the final time step is always saved.

- **Auto-scaling `max_iter`**: Poisson solver iterations now scale
  automatically with grid size: `max(1000, 20 * max(nx, ny))`. This
  ensures adequate iterations for larger grids (e.g., 128×128 now uses
  2560 iterations by default).

- **Suppressed initial transient warnings**: Poisson solver convergence
  warnings are now suppressed for the first 10 time steps, as the
  initial transient commonly has larger residuals that resolve quickly.

- **Optimized Jacobi solver memory**: Jacobi iterations now use
  pre-allocated buffer instead of cloning the matrix each iteration,
  reducing memory allocation overhead.

#### Documentation

- Updated Taylor-Green vortex formulas in all vignettes to match package
  sign convention.
- Clarified first-order accuracy (O(Δx) spatial from upwind, O(Δt)
  temporal from Euler).
- Added numerical diffusion quantification (~\|u\|Δx/2 from upwind
  scheme).
- Added comprehensive OpenMP parallelization documentation.
- Softened performance benchmark claims.

#### Tests

- Added Taylor-Green energy decay verification test.
- Added boundary conditions at t=0 verification test.
- Added IC/forcing consistency test.
- Added domain scaling test for initial conditions.

------------------------------------------------------------------------

## navieRstokes 0.1.0

### Initial release

#### Features

- Complete 2D incompressible Navier-Stokes solver using Chorin’s
  projection method

- Highly optimized Rcpp computational kernels achieving ~500x
  performance improvement:

  - Convection computation with upwind scheme
  - Diffusion (Laplacian) operator
  - Divergence and gradient operators
  - Optimized Jacobi Poisson equation solver with periodic/Neumann BC
    support
  - Vorticity calculation
  - Boundary condition enforcement

- Multiple boundary condition types:

  - Dirichlet (specified values)
  - Neumann (zero gradient)
  - Periodic (fully working with perfect mass conservation)

- Iterative pressure solvers:

  - Jacobi method (default, optimized for stability)
  - SOR (Successive Over-Relaxation) method available
  - Automatic adaptation for non-uniform grids
  - Support for periodic and Neumann boundary conditions in pressure
  - Fixed pressure pinning to remove singularity in pure Neumann BCs

- Comprehensive diagnostics:

  - Mass conservation error tracking
  - CFL number monitoring (diffusion and convection)
  - Pressure solver convergence tracking
  - Maximum velocity tracking

- Built-in visualization functions:

  - S3 methods for plot() and flow()
  - Velocity field plots with vector overlays
  - Vorticity contour plots
  - Pressure distribution plots
  - Animation support (MP4 and GIF) via av and gifski packages
  - A comprehensive Shiny dashboard simulating several scenarios.

- External forcing support for simulating obstacles and body forces

- Extensive documentation with four comprehensive vignettes:

  - Introduction to navieRstokes
  - Lid-driven cavity flow (classic benchmark)
  - Vortex decay simulation
  - Advanced examples (obstacle, channel flow, Taylor-Green,
    Kelvin-Helmholtz, etc.)

#### Performance

**v0.1.0 benchmarks** (optimized Jacobi solver):

- **64×64 grid, 1000 steps**: ~0.4 seconds
- **128×128 grid, 1000 steps**: ~2.5 seconds
- **256×256 grid, 1000 steps**: ~30 seconds

Overall speedup: ~500x compared to early development versions

#### Implementation details

- Second-order accurate spatial discretization
- Explicit time integration
- Upwind differencing for convection (first-order accurate, numerically
  stable)
- Central differences for diffusion, divergence, and gradient
  (second-order accurate)
- Automatic CFL condition checking with warnings
- S3 class system for clean interface
- Fixed Poisson solver singularity with pure Neumann BCs by pinning
  pressure at reference point

#### Solver improvements

- **Critical fix**: Resolved Poisson solver convergence issue
  - Root cause: Pure Neumann boundary conditions created singular system
  - Solution: Pin pressure at one point (`p[0,0] = 0`) to remove
    singularity
  - Switched default pressure solver from SOR to Jacobi (more stable)
  - Optimized default parameters: `max_iter = 1000`, `tolerance = 1e-3`
  - **Result**: ~500x speedup (64×64 grid: 221s → 0.4s for 1000 steps)

#### What works OK

- **Package compilation and installation**: All Rcpp files compile
  correctly, package installs without errors
- **Dirichlet boundary conditions**: Fast and stable with optimized
  Jacobi solver
  - Lid-driven cavity benchmark: 0.4s for 64×64 grid, 1000 steps
  - Excellent mass conservation with proper time stepping
- **Periodic boundary conditions**: Perfect mass conservation (error =
  0.00e+00)
  - Pressure Poisson solver with periodic BCs works correctly
  - Uniform flow test passes perfectly
  - Vortex decay simulations work flawlessly
- **CFL monitoring**: Both diffusion and convection CFL tracked
  correctly
- **All Rcpp optimizations**: 7 computational kernels fully functional

#### Known issues

- **SOR solver convergence**: SOR method has convergence problems -
  recommend using Jacobi (default)
- **Mass conservation error**: Can be relatively high (~6e-2) with CFL
  violations or improper time stepping
- **Periodic BCs**: Require divergence-free initial conditions for
  stability
- **CFL violations**: When dt is too large, simulations can become
  unstable
  - Solution: Package warns about CFL violations; users must respect
    these conditions

#### Test results

| Test | Boundary conditions | Result | Mass error | Performance |
|----|----|----|----|----|
| Uniform flow | Periodic | ✅ PASS | 0.00e+00 | OK |
| Vortex decay | Periodic | ✅ PASS | ~0.00e+00 | OK |
| Lid-driven cavity (Re=100) | Dirichlet | ✅ PASS | ~1e-3 | 0.4s (64×64, 1000 steps) |

#### Recommendations for users

1.  **Use default parameters** (Jacobi solver) for best performance and
    stability
2.  **For Dirichlet BCs**: Excellent stability and speed with proper
    time stepping
3.  **For periodic BCs**: Use package-provided divergence-free initial
    conditions with small amplitude
4.  **Respect CFL conditions**: The package warns you; listen to the
    warnings
5.  **Time step selection**: Follow diffusion CFL: `dt < dx²/(4*nu)` and
    convection CFL: `dt < dx/u_max`

#### Known limitations

- 2D only (no 3D implementation)
- Uniform Cartesian grid (no adaptive mesh refinement)
- First-order time integration
- Collocated grid (potential for pressure checkerboard modes)
- No turbulence modeling

#### Future development

Planned features for future releases:

- Higher-order time integration schemes (Runge-Kutta, Crank-Nicolson)
- Staggered grid option (MAC grid)
- 3D implementation
- Adaptive mesh refinement
- Multigrid Poisson solver for even better performance
- Turbulence modeling (LES, RANS)
- Heat transfer coupling
- Non-Newtonian fluid models
