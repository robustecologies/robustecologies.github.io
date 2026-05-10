# Compile a parallel basin computation function for user-defined dynamics

Generates and compiles a complete C++ source file containing the user's
dynamics function, an RK4 integrator, and an OpenMP-parallelized grid
computation loop. This provides performance matching the built-in
systems.

## Usage

``` r
compile_basin_function(
  cpp_dynamics,
  params = list(),
  dim = 2L,
  attractors,
  type = c("ode", "map"),
  verbose = FALSE
)
```

## Arguments

- cpp_dynamics:

  Character string containing C++ code for the dynamics. The code should
  compute derivatives and store them in a pre-allocated array. Available
  variables:

  - `state[i]` - Current state (read)

  - `deriv[i]` - Derivatives to compute (write)

  - `t` - Current time

  - `dim` - State dimension

  - Parameters by name (e.g., `damping`, `forcing`)

- params:

  Named list of parameters. These become C++ constants in the compiled
  code for maximum performance.

- dim:

  Integer. State space dimension.

- attractors:

  List of attractor specifications (from
  [`attractor_point()`](https://robustecologies.github.io/wadaR/reference/attractor_point.md),
  etc.)

- type:

  Character. System type: "ode" or "map".

- verbose:

  Logical. Show compilation output.

## Value

A function `f(x_grid, y_grid, t_max, dt, n_cores)` that computes basins
in parallel and returns an integer matrix.

## Details

The generated C++ code includes:

- User dynamics function (inline for performance)

- RK4 integrator with attractor convergence checking

- OpenMP parallel loop over grid points

- Periodic interrupt checking for Esc cancellation

**Example dynamics code for Duffing oscillator:**


    deriv[0] = state[1];
    deriv[1] = -delta * state[1] - alpha * state[0]
               - beta * state[0] * state[0] * state[0]
               + gamma_f * cos(omega * t);

## References

Butcher, J. C. (2008). *Numerical Methods for Ordinary Differential
Equations* (2nd ed.). Wiley. ISBN 978-0-470-72335-7.

Dagum, L., & Menon, R. (1998). OpenMP: an industry-standard API for
shared-memory programming. *IEEE Computational Science and Engineering*,
5(1), 46-55. [doi:10.1109/99.660313](https://doi.org/10.1109/99.660313)

Eddelbuettel, D., & Francois, R. (2011). Rcpp: Seamless R and C++
integration. *Journal of Statistical Software*, 40(8), 1-18.
[doi:10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08)

## See also

[`compiled_system`](https://robustecologies.github.io/wadaR/reference/compiled_system.md),
[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md),
[`attractor_point`](https://robustecologies.github.io/wadaR/reference/attractor_point.md),
[`attractor_cycle`](https://robustecologies.github.io/wadaR/reference/attractor_cycle.md),
[`is.compiled_system`](https://robustecologies.github.io/wadaR/reference/is.compiled_system.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ============================================================
# Example 1: Duffing oscillator (double-well potential)
# ============================================================
# The Duffing equation: x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t)

compute_duffing <- compile_basin_function(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -delta * state[1] - alpha * state[0]
                   - beta * state[0] * state[0] * state[0]
                   + gamma_f * cos(omega * t);
    ',
    params = list(delta = 0.3, alpha = -1, beta = 1,
                  gamma_f = 0.37, omega = 1.2),
    dim = 2,
    attractors = list(
        attractor_point(c(1, 0), 0.3, "Right well"),
        attractor_point(c(-1, 0), 0.3, "Left well"),
        attractor_point(c(0, 0), 0.3, "Origin")
    ),
    verbose = TRUE
)

# Compute basins with full OpenMP parallelization
basins <- compute_duffing(
    x_grid = seq(-2, 2, length.out = 300),
    y_grid = seq(-2, 2, length.out = 300),
    t_max = 100,
    dt = 0.01,
    n_cores = 0  # Auto-detect
)

# Result is an integer matrix of basin assignments
image(basins)

# ============================================================
# Example 2: Van der Pol oscillator with forcing
# ============================================================
# x'' - mu*(1-x^2)*x' + x = A*cos(omega*t)

compute_vdp <- compile_basin_function(
    cpp_dynamics = '
        double x = state[0];
        double v = state[1];
        deriv[0] = v;
        deriv[1] = mu * (1.0 - x*x) * v - x + A * cos(omega * t);
    ',
    params = list(mu = 1.0, A = 1.2, omega = 1.0),
    dim = 2,
    attractors = list(
        attractor_cycle(c(2, 0), 0.5, label = "Large orbit"),
        attractor_cycle(c(-2, 0), 0.5, label = "Small orbit"),
        attractor_point(c(0, 0), 0.3, label = "Origin")
    )
)

# ============================================================
# Example 3: Henon map (discrete iteration)
# ============================================================
# x_{n+1} = 1 - a*x_n^2 + y_n
# y_{n+1} = b*x_n

compute_henon <- compile_basin_function(
    cpp_dynamics = '
        next_state[0] = 1.0 - a * state[0] * state[0] + state[1];
        next_state[1] = b * state[0];
    ',
    params = list(a = 1.4, b = 0.3),
    dim = 2,
    attractors = list(
        attractor_point(c(0.63, 0.19), 0.15, "Fixed point 1"),
        attractor_point(c(-1.13, -0.34), 0.15, "Fixed point 2"),
        attractor_point(c(0, 0), 5.0, "Escape")
    ),
    type = "map"
)

henon_basins <- compute_henon(
    x_grid = seq(-1.5, 1.5, length.out = 400),
    y_grid = seq(-0.5, 0.5, length.out = 400),
    t_max = 500,  # Max iterations
    dt = 1,       # Ignored for maps
    n_cores = 0
)
} # }
```
