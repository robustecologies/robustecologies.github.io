# Changelog

## wadaR 1.1.0

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

------------------------------------------------------------------------

## wadaR 1.0.0

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

### References

1.  Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
    Nonlinear Phenomena*, 51(1-3), 213-225.
    <doi:10.1016/0167-2789(91)90234-Z>

2.  Nusse, H. E., & Yorke, J. A. (1996). Wada basin boundaries and basin
    cells. *Physica D: Nonlinear Phenomena*, 90(3), 242-261.
    <doi:10.1016/0167-2789(95)00249-9>

3.  Daza, A., Wagemakers, A., Sanjuan, M. A. F., & Yorke, J. A. (2015).
    Testing for basins of Wada. *Scientific Reports*, 5, 16579.
    <doi:10.1038/srep16579>

4.  Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2018). Ascertaining
    when a basin is Wada: The merging method. *Scientific Reports*,
    8, 9954. <doi:10.1038/s41598-018-28119-0>

5.  Wagemakers, A., Daza, A., & Sanjuan, M. A. F. (2020). The
    saddle-straddle method to test for Wada basins. *Communications in
    Nonlinear Science and Numerical Simulation*, 84, 105167.
    <doi:10.1016/j.cnsns.2020.105167>

6.  Wagemakers, A., Daza, A., & Sanjuan, M. A. F. (2021). How to detect
    Wada basins. *Discrete and Continuous Dynamical Systems - B*, 26(1),
    717-739. <doi:10.3934/dcdsb.2020330>
