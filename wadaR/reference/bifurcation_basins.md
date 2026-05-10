# Bifurcation analysis of basins of attraction

Performs parameter sweeps to analyze how basins of attraction change as
a system parameter varies. Uses C++/OpenMP parallelization for high
performance with interrupt support.

## Usage

``` r
bifurcation_basins(
  cpp_dynamics,
  params,
  sweep_param,
  sweep_values,
  dim = 2,
  attractors,
  type = c("ode", "map"),
  x_range,
  y_range,
  resolution = 100,
  t_max = 100,
  dt = 0.01,
  n_cores = NULL,
  verbose = TRUE,
  compute_entropy = TRUE,
  compute_wada = FALSE
)
```

## Arguments

- cpp_dynamics:

  Character string containing C++ code for the dynamics. Same format as
  [`compile_basin_function`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md).

- params:

  Named list of fixed parameters.

- sweep_param:

  Character. Name of the parameter to sweep.

- sweep_values:

  Numeric vector. Values of the sweep parameter.

- dim:

  Integer. State space dimension (default 2).

- attractors:

  List of attractor specifications or a function that takes the sweep
  parameter value and returns a list of attractors. This allows
  attractors to move with the parameter.

- type:

  Character. System type: "ode" or "map" (default "ode").

- x_range:

  Numeric vector of length 2. Range of x values.

- y_range:

  Numeric vector of length 2. Range of y values.

- resolution:

  Integer. Grid resolution (default 100).

- t_max:

  Numeric. Maximum integration time (default 100).

- dt:

  Numeric. Integration time step (default 0.01).

- n_cores:

  Integer. Number of CPU cores to use for parallel computation. Default
  is `NULL`, which uses `parallel::detectCores(logical = FALSE) - 1`.

- verbose:

  Logical. Show progress (default TRUE).

- compute_entropy:

  Logical. Compute basin entropy at each step (default TRUE).

- compute_wada:

  Logical. Compute Wada statistic at each step (default FALSE).

## Value

A list of class `bifurcation_result` containing:

- param_name:

  Name of the swept parameter.

- param_values:

  Vector of parameter values.

- basins:

  List of basin matrices, one per parameter value.

- entropy:

  Numeric vector of basin entropy values (if computed).

- boundary_fraction:

  Numeric vector of boundary fractions.

- n_basins:

  Integer vector of number of distinct basins detected.

- wada_statistic:

  Numeric vector of W statistics (if computed).

- x_grid:

  The x grid used.

- y_grid:

  The y grid used.

- attractors:

  The attractor specifications used.

## Details

Bifurcation analysis reveals how the structure of basins changes with
parameter variation. Key phenomena that can be detected include basin
erosion (gradual reduction in basin area), basin metamorphosis
(qualitative changes in basin shape), crisis events (sudden changes in
attractor/basin structure), and Wada transitions (emergence or loss of
the Wada property).

The computation is fully parallelized using OpenMP. Each parameter value
triggers a complete basin computation, and progress can be interrupted
with Esc.

## References

- Aguirre, J., Viana, R. L., & Sanjuan, M. A. F. (2009). Fractal
  structures in nonlinear dynamics. *Reviews of Modern Physics*, 81(1),
  333-386.
  [doi:10.1103/RevModPhys.81.333](https://doi.org/10.1103/RevModPhys.81.333)

- Grebogi, C., Ott, E., & Yorke, J. A. (1983). Crises, sudden changes in
  chaotic attractors, and transient chaos. *Physica D*, 7(1-3), 181-200.
  [doi:10.1016/0167-2789(83)90126-4](https://doi.org/10.1016/0167-2789%2883%2990126-4)

## See also

[`compile_basin_function`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md),
[`basin_entropy`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md),
[`plot.bifurcation_result`](https://robustecologies.github.io/wadaR/reference/plot.bifurcation_result.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ============================================================ #
# Example 1: Duffing oscillator - forcing amplitude sweep      #
# ============================================================ #
# Study how basins change as forcing increases

bif_duffing <- bifurcation_basins(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -delta * state[1] - alpha * state[0]
                   - beta * pow(state[0], 3)
                   + gamma_f * cos(omega * t);
    ',
    params = list(delta = 0.3, alpha = -1, beta = 1, omega = 1.2),
    sweep_param = "gamma_f",
    sweep_values = seq(0.1, 0.5, length.out = 20),
    attractors = list(
        attractor_point(c(1, 0), 0.3, "Right well"),
        attractor_point(c(-1, 0), 0.3, "Left well"),
        attractor_exit(0, "Escape")
    ),
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    resolution = 150,
    t_max = 100
)

# Plot bifurcation diagram
plot(bif_duffing)

# Animate basin evolution
animate(bif_duffing, filename = "duffing_bifurcation.gif")

# ============================================================ #
# Example 2: Henon map - parameter b sweep                     #
# ============================================================ #

bif_henon <- bifurcation_basins(
    cpp_dynamics = '
        next_state[0] = 1.0 - a * state[0] * state[0] + state[1];
        next_state[1] = b * state[0];
    ',
    params = list(a = 1.4),
    sweep_param = "b",
    sweep_values = seq(0.1, 0.4, length.out = 15),
    type = "map",
    attractors = list(
        attractor_point(c(0.63, 0.19), 0.2, "Attractor 1"),
        attractor_point(c(-1.13, -0.34), 0.2, "Attractor 2"),
        attractor_exit(0, "Escape")
    ),
    x_range = c(-2, 2),
    y_range = c(-1, 1),
    resolution = 200,
    t_max = 1000  # iterations for maps
)

plot(bif_henon, type = "entropy")
} # }
```
