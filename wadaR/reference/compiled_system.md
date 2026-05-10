# Create a compiled dynamical system with parallel basin computation

Creates a `compiled_system` object with pre-compiled C++/OpenMP parallel
basin computation. This is the recommended way to define custom systems
for high-performance Wada analysis.

## Usage

``` r
compiled_system(
  cpp_dynamics = NULL,
  model = NULL,
  attractors,
  dim = 2L,
  type = c("ode", "map"),
  params = list(),
  name = "Compiled system",
  description = NULL,
  verbose = FALSE
)
```

## Arguments

- cpp_dynamics:

  Character string containing C++ code for the dynamics. For ODEs,
  compute derivatives: `deriv[0] = ...; deriv[1] = ...;` For maps,
  compute next state: `next_state[0] = ...; next_state[1] = ...;`
  Default `NULL`; either `cpp_dynamics` or `model` must be supplied, not
  both.

- model:

  Optional
  [`janos::model_spec`](https://robustecologies.github.io/janos/reference/model_spec.html)
  object built from a formula list. When supplied, the per-state
  right-hand-side expressions are extracted via
  [`janos::model_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/model_spec_rhs_cpp.html)
  and rewritten into wadaR's C++ template (state references as
  `state[i]`, parameter references as bare names matched by the
  `const double` declarations that
  [`compile_basin_function()`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md)
  emits from the `params` list). Only deterministic ODE or discrete-map
  model_spec objects are accepted; stochastic, delayed, jump, spatial,
  switched and Markov-chain variants raise an error. When `model` is
  supplied, `dim`, `type` and `params` default to the values stored on
  the model_spec, but any explicit argument the user passes wins.
  Mutually exclusive with `cpp_dynamics`. Requires the `janos` package
  (in `Suggests`). Default: `NULL`.

- attractors:

  List of attractor specifications created with
  [`attractor_point()`](https://robustecologies.github.io/wadaR/reference/attractor_point.md),
  [`attractor_cycle()`](https://robustecologies.github.io/wadaR/reference/attractor_cycle.md),
  etc.

- dim:

  Integer. State space dimension (default 2).

- type:

  Character. System type: "ode" (default) or "map".

- params:

  Named list of scalar numeric parameters. These become compile-time
  constants in C++ for maximum performance.

- name:

  Character. Human-readable system name.

- description:

  Character. System description.

- verbose:

  Logical. Show compilation output.

## Value

A `compiled_system` object with a pre-compiled parallel basin
computation function stored in `$compiled_basin_func`.

## Details

This function generates and compiles optimized C++ code for your
dynamical system. The compiled code includes:

- Inline dynamics function (no function call overhead)

- Inline RK4 integrator

- OpenMP parallelization matching built-in systems

- Esc interrupt support for cancellation

**Performance**: Compiled systems run at the same speed as the built-in
pendulum and Henon-Heiles implementations.

**ODE dynamics code format:**


    // Available: state[i], t, deriv[i], dim, and your parameters
    deriv[0] = state[1];                    // dx/dt = v
    deriv[1] = -damping * state[1] - ...;   // dv/dt = ...

**Map dynamics code format:**


    // Available: state[i], iter, next_state[i], dim, and your parameters
    next_state[0] = 1 - a * state[0] * state[0] + state[1];
    next_state[1] = b * state[0];

**Alternative input via janos::model_spec.** A user who has already
built a `model_spec` for trajectory simulation in the companion package
`janos` can pass it through the `model` argument instead of writing the
dynamics in raw C++. The two routes are numerically equivalent: wadaR
consumes the same parsed expressions that janos's RK4/RK45 templates
would inline, rewrites them into its basin-kernel identifier convention,
and compiles the same parallel OpenMP grid sweep. Only formula-based
deterministic ODE and discrete-map model_spec objects are accepted;
stochastic, delayed, jump, spatial, switched and Markov-chain families
do not have an inlining contract compatible with the basin kernel and
raise an error.

## References

Eddelbuettel, D., & Francois, R. (2011). Rcpp: Seamless R and C++
integration. *Journal of Statistical Software*, 40(8), 1-18.
[doi:10.18637/jss.v040.i08](https://doi.org/10.18637/jss.v040.i08)

Dagum, L., & Menon, R. (1998). OpenMP: an industry-standard API for
shared-memory programming. *IEEE Computational Science and Engineering*,
5(1), 46-55. [doi:10.1109/99.660313](https://doi.org/10.1109/99.660313)

Butcher, J. C. (2008). *Numerical Methods for Ordinary Differential
Equations* (2nd ed.). Wiley. ISBN 978-0-470-72335-7.

## See also

[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md),
[`compile_basin_function`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md),
[`attractor_point`](https://robustecologies.github.io/wadaR/reference/attractor_point.md),
[`attractor_cycle`](https://robustecologies.github.io/wadaR/reference/attractor_cycle.md),
[`attractor_exit`](https://robustecologies.github.io/wadaR/reference/attractor_exit.md),
[`attractor_outcome`](https://robustecologies.github.io/wadaR/reference/attractor_outcome.md),
[`is.compiled_system`](https://robustecologies.github.io/wadaR/reference/is.compiled_system.md),
[`print.compiled_system`](https://robustecologies.github.io/wadaR/reference/print.compiled_system.md),
[`summary.compiled_system`](https://robustecologies.github.io/wadaR/reference/summary.compiled_system.md),
[`plot.compiled_system`](https://robustecologies.github.io/wadaR/reference/plot.compiled_system.md);
the formula-based interoperability route uses
[`janos::model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.html)
and
[`janos::model_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/model_spec_rhs_cpp.html).

## Examples

``` r
if (FALSE) { # \dontrun{
# ============================================================
# Example 1: Duffing oscillator - complete workflow
# ============================================================
# The Duffing oscillator with double-well potential exhibits
# chaotic behavior and fractal basin boundaries.
#
# Equation: x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t)

duffing <- compiled_system(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -delta * state[1] - alpha * state[0]
                   - beta * state[0] * state[0] * state[0]
                   + gamma_f * cos(omega * t);
    ',
    attractors = list(
        attractor_point(c(1, 0), 0.3, "Right well"),
        attractor_point(c(-1, 0), 0.3, "Left well"),
        attractor_point(c(0, 0), 0.3, "Origin")
    ),
    params = list(delta = 0.3, alpha = -1, beta = 1,
                  gamma_f = 0.37, omega = 1.2),
    name = "Duffing oscillator",
    description = "Forced Duffing with double-well potential"
)

# Print system information
print(duffing)

# Check if system is compiled
is.compiled_system(duffing)  # TRUE

# Compute basins with OpenMP parallelization
basins <- compute_basins(
    duffing,
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    resolution = 300,
    t_max = 100,
    dt = 0.01,
    verbose = TRUE
)

# Visualize basins
plot(basins, title = "Duffing oscillator basins")

# Compute basin entropy
entropy <- basin_entropy(basins, box_size = 2)
print(entropy)
plot(entropy)

# Test for Wada property
grid_result <- wada_grid_method(basins, verbose = TRUE)
print(grid_result)

merge_result <- wada_merging_method(basins, verbose = TRUE)
print(merge_result)

# ============================================================
# Example 2: Henon map - discrete dynamics
# ============================================================
# Classic area-preserving map with chaotic dynamics.
# x_{n+1} = 1 - a*x_n^2 + y_n
# y_{n+1} = b*x_n

henon <- compiled_system(
    cpp_dynamics = '
        next_state[0] = 1.0 - a * state[0] * state[0] + state[1];
        next_state[1] = b * state[0];
    ',
    attractors = list(
        attractor_point(c(0.63, 0.19), 0.15, "Fixed point 1"),
        attractor_point(c(-1.13, -0.34), 0.15, "Fixed point 2"),
        attractor_point(c(0, 0), 5.0, "Escape")
    ),
    type = "map",
    params = list(a = 1.4, b = 0.3),
    name = "Henon map"
)

# Compute basins (max iterations for map)
basins_henon <- compute_basins(
    henon,
    x_range = c(-1.5, 1.5),
    y_range = c(-0.5, 0.5),
    resolution = 400,
    t_max = 500  # Max iterations
)

plot(basins_henon, title = "Henon map basins")

# ============================================================
# Example 3: Van der Pol oscillator with forcing
# ============================================================
# Self-sustained oscillator with limit cycle behavior.
# x'' - mu*(1-x^2)*x' + x = A*cos(omega*t)

vdp <- compiled_system(
    cpp_dynamics = '
        double x = state[0];
        double v = state[1];
        deriv[0] = v;
        deriv[1] = mu * (1.0 - x*x) * v - x + A * cos(omega * t);
    ',
    attractors = list(
        attractor_cycle(c(2, 0), 0.5, label = "Large orbit"),
        attractor_cycle(c(-2, 0), 0.5, label = "Small orbit"),
        attractor_point(c(0, 0), 0.3, label = "Origin")
    ),
    params = list(mu = 1.0, A = 1.2, omega = 1.0),
    name = "Van der Pol oscillator"
)

print(vdp)

# ============================================================
# Example 4: Same Duffing system through janos::model_spec()
# ============================================================
# When janos is installed, dynamics can be specified through
# formulas instead of raw C++; the resulting compiled_system is
# numerically identical to the raw-C++ route.

if (requireNamespace("janos", quietly = TRUE)) {
    duffing_via_janos <- compiled_system(
        model = janos::model_spec(
            rhs = list(
                x ~ y,
                y ~ -delta * y - alpha * x - beta * x * x * x
                    + gamma_f * cos(omega * t)
            ),
            state_names = c("x", "y"),
            parms = list(delta = 0.3, alpha = -1, beta = 1,
                         gamma_f = 0.37, omega = 1.2)
        ),
        attractors = list(
            attractor_point(c( 1, 0), 0.3, "Right well"),
            attractor_point(c(-1, 0), 0.3, "Left well"),
            attractor_point(c( 0, 0), 0.3, "Origin")
        ),
        name = "Duffing oscillator (via janos::model_spec)"
    )
    print(duffing_via_janos)
}
} # }
```
