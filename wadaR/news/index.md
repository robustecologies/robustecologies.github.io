# Changelog

## wadaR 0.3.3

### API changes

- All twelve exported `plot.*` S3 methods gain a `caption` argument
  (default `TRUE`). When `TRUE`, the caption is a tightened one-line
  string of at most 80 characters carrying the source function name and
  primary citation; when `FALSE`, the caption is suppressed entirely.
  This eliminates caption overlap and truncation under faceted,
  patchwork-composed or small-device rendering. Run-specific numerics
  that previously lived in the caption (max d_H, relative diff,
  slim-boundary point count) now appear in the subtitle. The
  `plot.wada_analysis` and `plot.bifurcation_result` patchwork branches
  automatically suppress per-panel captions and render a single caption
  on the outer annotation. The 3D `plot.basin_result_3d` accepts the
  argument for API uniformity; it has effect only on the 2D-slice
  branch, since plotly widgets carry no caption layer.

  

## wadaR 0.3.2

### New features

- [`palettes()`](https://robustecologies.github.io/wadaR/reference/palettes.md)
  and
  [`palette_ramp()`](https://robustecologies.github.io/wadaR/reference/palette_ramp.md)
  now delegate to the koloRo package when installed, exposing the full
  koloRo catalogue (282 palettes across scientific, colorblind-safe,
  alhambra, chameleons, natural, cultural, artistic and seasonal
  categories) directly from wadaR. koloRo is declared in `Suggests`, so
  wadaR remains installable without it; the colorblind-safe `okabe_ito`
  palette is retained as a built-in fallback so the bundled
  `wada_basins`, `wada_straddle_method` and `wada_merging_method` plots
  continue to render unchanged. Any other palette name requested without
  koloRo installed raises an explicit error with the install
  instruction.

- [`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
  accepts an optional `model` argument: a
  [`janos::model_spec`](https://robustecologies.github.io/janos/reference/model_spec.html)
  object built from formula-based dynamics is translated transparently
  into wadaR’s C++ basin kernel via the new
  [`janos::model_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/model_spec_rhs_cpp.html)
  helper (Suggests). This avoids restating the equations in raw C++ when
  the system has already been specified for
  [`janos::dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.html).
  State references in the user formulas are rewritten to wadaR’s
  `state[i]` convention and parameter references to bare names matched
  by the `const double` declarations of the basin template; the
  resulting compiled object is numerically identical to the raw-C++
  route. Only deterministic ODE and discrete-map model_spec objects are
  supported; stochastic, delayed, jump, spatial, switched and
  Markov-chain variants raise an explicit error.

- [`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
  argument `cpp_dynamics` is now optional (default `NULL`). One of
  `cpp_dynamics` or `model` must be supplied; supplying both raises an
  error. Existing scripts that pass `cpp_dynamics` positionally or by
  name continue to work unchanged.

## wadaR 0.3.1

### API hygiene and CRAN policy compliance

- All [`stop()`](https://rdrr.io/r/base/stop.html) calls in exported
  functions now pass `call. = FALSE`, so error messages no longer leak
  the internal call stack into the user-facing report.

- Plot S3 methods on every ggplot2-backed class now carry a synthetic
  technical subtitle and a grey caption that names the source function
  and the relevant numeric convention. The plotly methods on
  `basin_result_3d` retain their widget-native annotations.

- [`detect_wada()`](https://robustecologies.github.io/wadaR/reference/detect_wada.md)
  now documents its dispatch contract explicitly, expands the
  `@references` block to cite all three method papers (Daza et al. 2015,
  Daza et al. 2018, Wagemakers et al. 2020) plus the foundational
  Kennedy and Yorke (1991), and bidirectionally cross-references the
  `wada_analysis` S3 trio.

- [`compile_basin_function()`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md)
  and
  [`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
  gain `@references` to Butcher (2008), Dagum and Menon (1998), and
  Eddelbuettel and Francois (2011), making the underlying
  RK4/OpenMP/Rcpp foundations citable from the help pages.

- [`is.compiled_system()`](https://robustecologies.github.io/wadaR/reference/is.compiled_system.md)
  is upgraded from a one-line predicate to a full roxygen entry with
  `@details`, runnable example and bidirectional `@seealso`.

- The package now declares a minimum R version of 4.1.0, reflecting the
  use of native pipe syntax in the 3D visualization layer. R 4.1.0 has
  been the universally available baseline since 2021.

### Tests

- New `tests/testthat/setup.R` pins `OMP_NUM_THREADS=1` for reproducible
  parallel-Rcpp test results on multi-core hosts.

- New regression tests for `bifurcation_basins`, `extract_basin`,
  `compute_basins_3d`, `slice_3d_basins`, `plot_3d_basins`,
  `detect_wada`, `palettes`, `palette_ramp`, the `animate` generic and
  `shinywadaR`. Long-running tests use `skip_on_cran()` and
  Suggests-dependent tests use `skip_if_not_installed()`.

- The `hausdorff_distance` symmetry test now seeds `set.seed(20260425)`
  to remove a residual stochastic dependency.

  

## wadaR 0.3.0

### Validation guards and canonical regression tests

- `wada_grid_method`, `wada_merging_method` and `wada_straddle_method`
  now error early with the message
  `"Wada detection requires at least 3 attractors; got N."` when called
  on basins with fewer than three attractors. Previously the methods
  would proceed and either return a numerically meaningless result or
  fail with a cryptic Rcpp error.

- New canonical-validation test suite
  (`tests/testthat/test-validation-canonical.R`, `skip_on_cran()`) pins
  the Wada classification on three published reference systems: the
  Newton fractal for `z^3 - 1` (Wada at 128x128 grid), the forced damped
  pendulum at `F = 1.66` (Wada) and at `F = 2.5` (non-Wada control),
  plus the Huisman-Weissing competition systems (5sp -\> 2 attractors,
  8sp -\> 3 attractors). Failure of these tests indicates an algorithmic
  regression and aborts release.

- `wada_merging_method`’s empirical thresholds (`relative_diff < 2`,
  `relative_max < 0.05`) are now documented explicitly in `@details` as
  defaults calibrated against the canonical benchmarks. Daza et
  al. (2018) propose adaptive thresholds tied to grid resolution; the
  constants are the package contract for the published benchmarks and
  should be re-validated when changing resolution or system class.

  

## wadaR 0.2.0

### New features

#### User-defined dynamical systems

- **[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)**:
  Create custom dynamical systems with C++/OpenMP compilation for
  maximum performance
  - Accepts dynamics as C++ code strings
  - Automatic OpenMP parallelization matching built-in systems
  - Full Esc interrupt support during computation
- **[`attractor_point()`](https://robustecologies.github.io/wadaR/reference/attractor_point.md)**:
  Helper function to define point attractors (fixed points)
- **[`attractor_cycle()`](https://robustecologies.github.io/wadaR/reference/attractor_cycle.md)**:
  Helper function to define limit cycle attractors
- **[`attractor_exit()`](https://robustecologies.github.io/wadaR/reference/attractor_exit.md)**:
  Helper function to define escape exits (for escape systems)
- **[`attractor_outcome()`](https://robustecologies.github.io/wadaR/reference/attractor_outcome.md)**:
  Helper function to define discrete outcomes (for competition systems)
- **[`is.compiled_system()`](https://robustecologies.github.io/wadaR/reference/is.compiled_system.md)**:
  Check if a system is compiled
- Users can now define arbitrary dynamical systems that work with all
  wadaR methods:
  [`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md),
  [`basin_entropy()`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md),
  [`wada_grid_method()`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md),
  [`wada_merging_method()`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md),
  etc.

#### Bifurcation analysis

- **[`bifurcation_basins()`](https://robustecologies.github.io/wadaR/reference/bifurcation_basins.md)**:
  Compute basins of attraction across a range of parameter values with
  C++/OpenMP parallelization
  - Supports any parameter in the dynamical system
  - Automatically computes basin entropy for each parameter value
  - Automatic boundary detection and Wada property estimation
  - Full user interrupt support (Esc to abort)
- **[`plot.bifurcation_result()`](https://robustecologies.github.io/wadaR/reference/plot.bifurcation_result.md)**:
  S3 plot method with multiple visualization modes
  - `mode = "summary"`: Basin entropy and boundary fraction
    vs. parameter
  - `mode = "basins"`: Faceted grid of basin plots across parameter
    values
  - `mode = "entropy"`: Local entropy maps for each parameter value
- **[`animate.bifurcation_result()`](https://robustecologies.github.io/wadaR/reference/animate.bifurcation_result.md)**:
  Create animated GIFs or HTML showing basin evolution across parameter
  values
- **[`extract_basin()`](https://robustecologies.github.io/wadaR/reference/extract_basin.md)**:
  Extract a single basin result at a specific parameter value as a
  `wada_basins` object
- **[`print.bifurcation_result()`](https://robustecologies.github.io/wadaR/reference/print.bifurcation_result.md)**:
  S3 print method with summary statistics

#### 3D visualization

- **[`compute_basins_3d()`](https://robustecologies.github.io/wadaR/reference/compute_basins_3d.md)**:
  Compute basins of attraction in 3D parameter space with C++/OpenMP
  parallelization
  - Supports any 3D projection of higher-dimensional systems
  - Full user interrupt support (Esc to abort)
- **[`plot_3d_basins()`](https://robustecologies.github.io/wadaR/reference/plot_3d_basins.md)**:
  Interactive 3D visualization using plotly
  - `mode = "scatter"`: Point cloud visualization
  - `mode = "isosurface"`: Surface rendering of basin boundaries
  - `mode = "volume"`: Volume rendering with opacity control
- **[`slice_3d_basins()`](https://robustecologies.github.io/wadaR/reference/slice_3d_basins.md)**:
  Extract 2D slices from 3D basin results
  - Slices along any axis (x, y, or z)
  - Returns `wada_basins` object compatible with all wadaR methods
- **[`animate_3d_rotation()`](https://robustecologies.github.io/wadaR/reference/animate_3d_rotation.md)**:
  Create animations of 3D basins
  - GIF output: Z-slice animation using gganimate
  - HTML output: Interactive 3D rotation using plotly
- **[`print.basin_result_3d()`](https://robustecologies.github.io/wadaR/reference/print.basin_result_3d.md)**:
  S3 print method with summary statistics

#### Runtime C++ compilation

- **[`compile_basin_function()`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md)**:
  Compile user-defined dynamics to C++ at runtime for maximum
  performance
- Uses
  [`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) for
  just-in-time compilation
- Automatic OpenMP parallelization of compiled code

#### New vignette

- **`vignettes/e.custom-systems.Rmd`**: Comprehensive guide to
  user-defined dynamical systems
  - Defining custom ODEs (Duffing oscillator example)
  - Bifurcation analysis with parameter sweeps
  - 3D visualization of higher-dimensional systems (Lorenz attractor)
  - Performance considerations and recommendations

### Improvements

- Enhanced verbose output with `get_colored_symbol()` for consistent
  terminal styling
- Improved computational complexity estimates displayed during
  computation
- All new functions support user interrupts (Esc to abort)

### Dependencies

- Added `gganimate` to Suggests for animated bifurcation plots
- Added `htmlwidgets` to Suggests for saving plotly animations

------------------------------------------------------------------------

## wadaR 0.1.0

### New features

#### Multispecies competition model

- **[`multispecies_competition()`](https://robustecologies.github.io/wadaR/reference/multispecies_competition.md)**:
  Factory function to create Huisman-Weissing multispecies competition
  systems with 5-species or 8-species scenarios
- **[`compute_competition_basins()`](https://robustecologies.github.io/wadaR/reference/compute_competition_basins.md)**:
  Compute basins of attraction for competition outcomes with fractal
  boundaries
- **[`simulate_competition()`](https://robustecologies.github.io/wadaR/reference/simulate_competition.md)**:
  Integrate competition dynamics and return time series

#### User-controlled parallelization

- **`n_cores` parameter**: All parallel functions now accept an
  `n_cores` parameter to control the number of CPU cores used
  - [`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md):
    Basin computation for pendulum and Henon-Heiles systems
  - [`compute_newton_basins()`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md):
    Newton fractal computation
  - [`compute_competition_basins()`](https://robustecologies.github.io/wadaR/reference/compute_competition_basins.md):
    Competition basin computation
  - All internal parallel functions (classify_balls, find_boundary,
    hausdorff_matrix, etc.)
- **Default**: `n_cores = NULL` uses
  `parallel::detectCores(logical = FALSE) - 1` (leaving one core free
  for system tasks)

#### User interrupt support (Esc to abort)

- **Automatic abort capability**: All computationally intensive parallel
  functions now check for user interrupts periodically
- Press **Esc** during computation to gracefully abort long-running
  simulations
- Uses `Rcpp::checkUserInterrupt()` with chunked processing to allow
  interrupt checking between parallel chunks

#### Basin entropy

- **[`basin_entropy()`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md)**:
  Compute basin entropy to quantify unpredictability in dynamical
  systems

#### New vignettes

- **`vignettes/competition-basins.Rmd`**: Comprehensive guide to fractal
  basins in multispecies competition, including ecological
  interpretation

### Improvements

- Enhanced verbose output showing number of cores being used and abort
  instructions
- Improved computational complexity estimates displayed during
  computation

### References

1.  Huisman, J., & Weissing, F. J. (2001). Fundamental unpredictability
    in multispecies competition. *American Naturalist*, 157(5), 488-494.
    <doi:10.1086/319929>

------------------------------------------------------------------------

## wadaR 0.0.9

### Major features

#### Wada basin detection methods

This package implements three complementary computational methods for
testing the Wada property in basins of attraction:

- **Grid method** (Daza et al., 2015): Tests whether a third color can
  always be found between any two colors at successively finer
  resolutions. Uses parallel Rcpp/OpenMP for high performance.

- **Merging method** (Daza et al., 2018): Exploits the property that
  Wada basins can be merged without changing their common boundary.
  Compares slim boundaries using the Hausdorff distance.

- **Saddle-straddle method** (Wagemakers et al., 2020): Tracks the
  chaotic saddle on basin boundaries and verifies that there is only one
  saddle (implying a single Wada boundary).

#### Example dynamical systems

- **Forced damped pendulum**: The canonical example for Wada basins with
  three coexisting period-1 attractors
- **Henon-Heiles system**: Hamiltonian system with escape basins
  exhibiting Wada boundaries
- **Newton fractals**: Mathematically proven Wada basins from
  Newton-Raphson iteration

#### High-performance implementations

- **OpenMP parallelization**: All computationally intensive operations
  parallelized with OpenMP
- **Rcpp integration**: Core algorithms implemented in C++ for maximum
  performance
- **Thread-safe design**: STL-based implementations avoid R’s
  PROTECT/UNPROTECT issues in parallel code

#### Visualization

- **S3 plot methods**: Consistent ggplot2-based visualization for all
  result types
- **Custom color palettes**: Full control over basin coloring
- **Multiple plot modes**: Overlay, faceted, and comparison views for
  different analysis needs

### Function summary

#### Basin computation

- [`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md):
  Compute basins of attraction for ODE systems
- [`compute_newton_basins()`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md):
  Specialized function for Newton fractals
- [`get_boundary()`](https://robustecologies.github.io/wadaR/reference/get_boundary.md):
  Extract boundary points from basin matrix
- [`merge_basins()`](https://robustecologies.github.io/wadaR/reference/merge_basins.md):
  Create two-color basin maps

#### System generators

- [`forced_damped_pendulum()`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md):
  Create pendulum system object
- [`henon_heiles_system()`](https://robustecologies.github.io/wadaR/reference/henon_heiles_system.md):
  Create Henon-Heiles system object
- [`newton_fractal_system()`](https://robustecologies.github.io/wadaR/reference/newton_fractal_system.md):
  Create Newton iteration system object

#### Wada detection

- [`wada_grid_method()`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md):
  Grid refinement method
- [`wada_merging_method()`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md):
  Boundary merging method
- [`wada_straddle_method()`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md):
  Saddle-straddle method
- [`detect_wada()`](https://robustecologies.github.io/wadaR/reference/detect_wada.md):
  Unified interface for multiple methods

#### Utilities

- [`hausdorff_distance()`](https://robustecologies.github.io/wadaR/reference/hausdorff_distance.md):
  Compute Hausdorff distance between point sets
- [`plot.wada_basins()`](https://robustecologies.github.io/wadaR/reference/plot.wada_basins.md):
  S3 plot method for basins
- [`plot.wada_grid_result()`](https://robustecologies.github.io/wadaR/reference/plot.wada_grid_result.md):
  S3 plot method for grid results
- [`plot.wada_merging_result()`](https://robustecologies.github.io/wadaR/reference/plot.wada_merging_result.md):
  S3 plot method for merging results
- [`plot.wada_straddle_result()`](https://robustecologies.github.io/wadaR/reference/plot.wada_straddle_result.md):
  S3 plot method for straddle results

### Test coverage

- 121 passing tests covering all major functions
- Tests for basin computation, S3 methods, system generators, Wada
  methods, and Newton fractals

### Dependencies

#### Required

- R (\>= 4.0.0)
- Rcpp (\>= 1.0.0)
- ggplot2 (\>= 3.4.0)
- R6

#### System requirements

- GNU make
- OpenMP (for parallel computation)

#### Suggested

- RANN (for fast nearest neighbor queries)
- testthat (\>= 3.0.0)
- knitr
- rmarkdown
