# Compute basins of attraction for a 2D dynamical system

Integrates trajectories from a grid of initial conditions to determine
which attractor each point converges to. Uses high-performance parallel
Rcpp implementations with OpenMP for all supported systems.

## Usage

``` r
compute_basins(
  system_func,
  x_range,
  y_range,
  resolution = 500,
  attractors = NULL,
  t_max = 100,
  dt = 0.01,
  is_map = FALSE,
  max_iter = 1000,
  system_type = c("auto", "pendulum", "henon_heiles"),
  params = NULL,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- system_func:

  Function(state, t) defining the dynamical system, OR a system object
  returned by
  [`forced_damped_pendulum()`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md),
  [`henon_heiles_system()`](https://robustecologies.github.io/wadaR/reference/henon_heiles_system.md),
  or
  [`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md).
  If a system object is passed, parameters are automatically extracted.
  For user-defined `compiled_system` objects, a compiled C++/OpenMP
  integrator is used for full parallel performance.

- x_range:

  Numeric vector of length 2. Range of x coordinates.

- y_range:

  Numeric vector of length 2. Range of y coordinates.

- resolution:

  Integer. Number of grid points per dimension.

- attractors:

  List of attractor specifications. Each element should be a list with
  'center' (numeric vector) and 'radius' (numeric). Not required if
  system_func is a system object with attractors.

- t_max:

  Numeric. Maximum integration time (for ODEs).

- dt:

  Numeric. Time step for integration.

- is_map:

  Logical. TRUE if system_func is a discrete map.

- max_iter:

  Integer. Maximum iterations for maps.

- system_type:

  Character. Type of system: "auto" (default, detects automatically),
  "pendulum", or "henon_heiles". All use Rcpp parallel.

- params:

  List of system-specific parameters (e.g., damping, forcing). Not
  required if system_func is a system object.

- n_cores:

  Integer. Number of CPU cores to use for parallel computation. Default
  is `NULL`, which uses `parallel::detectCores(logical = FALSE) - 1`
  (leaving one core free for system tasks). Set to 1 to disable
  parallelization.

- verbose:

  Logical. Show progress.

## Value

A list of class "wada_basins" containing:

- basins:

  Matrix of basin assignments (integers 1 to N)

- x_grid:

  Vector of x coordinates

- y_grid:

  Vector of y coordinates

- n_attractors:

  Number of attractors

- unclassified:

  Number of points that did not converge

## Details

Basins of attraction partition the phase space of a dynamical system
into regions that converge to different attractors. For dissipative
systems with multiple coexisting attractors, the basin structure can be
highly complex, exhibiting fractal boundaries and, in special cases, the
Wada property.

The algorithm uses fourth-order Runge-Kutta (RK4) integration to evolve
trajectories from each grid point until they enter a neighborhood of an
attractor or reach maximum integration time. The implementation is fully
parallelized using OpenMP for optimal performance on multi-core systems.

**Computational complexity warning:**

The runtime scales as \\O(N^2 \times T/\Delta t)\\ where \\N\\ is the
resolution, \\T\\ is `t_max`, and \\\Delta t\\ is `dt`. For a 500x500
grid with default parameters (`t_max=100`, `dt=0.01`), this means
approximately **2.5 billion RK4 steps**, which can take several minutes
even with parallel execution. Recommendations:

- Use resolution 100-200 for exploration

- Increase to 400-500 only for final publication-quality results

- Newton fractals (via
  [`compute_newton_basins`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md))
  are much faster as they don't require ODE integration

For the forced damped pendulum, the equation of motion is:
\$\$\ddot{x} + \gamma \dot{x} + \sin(x) = F \cos(t)\$\$ where \\\gamma\\
is the damping coefficient and \\F\\ is the forcing amplitude. This
system exhibits Wada basins for certain parameter values, notably around
\\F \approx 1.66\\ with \\\gamma = 0.2\\.

For the Henon-Heiles system, the Hamiltonian is: \$\$H =
\frac{1}{2}(p_x^2 + p_y^2) + \frac{1}{2}(x^2 + y^2) + x^2 y -
\frac{y^3}{3}\$\$ Trajectories can escape through three channels at
120-degree intervals when the energy exceeds the critical value \\E_c =
1/6\\.

## Note

**User abort**: Press `Esc` during computation to abort the simulation.
The function checks for user interrupts periodically and will stop
gracefully.

## References

Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
Nonlinear Phenomena*, 51(1-3), 213-225.
[doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789%2891%2990234-Z)

Aguirre, J., Vallejo, J. C., & Sanjuan, M. A. F. (2001). Wada basins and
chaotic invariant sets in the Henon-Heiles system. *Physical Review E*,
64(6), 066208.
[doi:10.1103/PhysRevE.64.066208](https://doi.org/10.1103/PhysRevE.64.066208)

Feudel, U. (2008). Complex dynamics in multistable systems.
*International Journal of Bifurcation and Chaos*, 18(06), 1607-1626.
[doi:10.1142/S0218127408021233](https://doi.org/10.1142/S0218127408021233)

## See also

[`compiled_system`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
for creating custom user-defined systems,
[`forced_damped_pendulum`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md)
for creating pendulum systems,
[`henon_heiles_system`](https://robustecologies.github.io/wadaR/reference/henon_heiles_system.md)
for creating Henon-Heiles systems,
[`compute_newton_basins`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md)
for Newton fractal basins,
[`plot.wada_basins`](https://robustecologies.github.io/wadaR/reference/plot.wada_basins.md)
for visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Forced damped pendulum with Wada basins                    #
# ===================================================================== #
# The forcing amplitude F = 1.66 produces clear Wada basins
pendulum <- forced_damped_pendulum(forcing = 1.66, damping = 0.2)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 500)
plot(basins, title = "Forced damped pendulum (F = 1.66)")

# ===================================================================== #
# Example 2: Partially Wada case                                        #
# ===================================================================== #
# F = 1.71 produces partially Wada basins (mixed Wada/non-Wada)
pendulum_partial <- forced_damped_pendulum(forcing = 1.71)
basins_partial <- compute_basins(pendulum_partial, c(-pi, pi), c(-3, 3),
                                 resolution = 500)
plot(basins_partial, title = "Partially Wada (F = 1.71)")

# ===================================================================== #
# Example 3: Henon-Heiles escape basins                                 #
# ===================================================================== #
# Energy above critical (E_c = 1/6) allows escape through 3 channels
hh <- henon_heiles_system(energy = 0.2)
basins_hh <- compute_basins(hh, c(-0.4, 0.4), c(-0.4, 0.4),
                            resolution = 500, t_max = 200)
plot(basins_hh, title = "Henon-Heiles escape basins (E = 0.2)",
     colors = c("#E41A1C", "#377EB8", "#4DAF4A"))

# ===================================================================== #
# Example 4: Show boundary structure                                    #
# ===================================================================== #
plot(basins, show_boundary = TRUE, boundary_color = "white",
     title = "Basins with fractal boundary overlay")
} # }
```
