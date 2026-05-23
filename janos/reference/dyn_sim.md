# Simulate a dynamical system

Integrates a dynamical system defined by a
[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md)
using the specified solver. This is the primary simulation entry point:
the model describes the dynamics while the solver determines the
integration method. Supports both ODE integration and stochastic
simulation of CTMC models.

## Usage

``` r
dyn_sim(
  model,
  t_max = 50,
  solver = solver_rk4(),
  init = NULL,
  parms = NULL,
  discard_transient = 30,
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- model:

  A `system_spec` object defining the dynamical system.

- t_max:

  Simulation duration in model time units (default: 50).

- solver:

  A `solver_spec` object (default:
  [`solver_rk4()`](https://robustecologies.github.io/janos/reference/solver_rk4.md)).

- init:

  Optional named numeric vector of initial conditions overriding the
  model defaults.

- parms:

  Optional named list of parameters overriding the model defaults.

- discard_transient:

  Time to discard as transient dynamics (default: 30).

- seed:

  Random seed for reproducible stochastic simulations (default: 42).

- verbose:

  Logical; print progress messages (default: TRUE).

## Value

An S3 object of class `dyn_sim` containing:

- trajectory:

  Full trajectory data frame with columns: time, state variables, and
  derived quantities

- attractor:

  Trajectory after discarding transients

- model:

  The `system_spec` used

- solver:

  The `solver_spec` used

- parameters:

  Named list of all parameters used (model + simulation)

- initial_conditions:

  Named vector of initial state

- metadata:

  Integration diagnostics: method, timing, step counts. For SSA solvers,
  includes `n_reactions` and `reaction_counts`.

## Details

The function dispatches to compiled C++ backends when the model provides
a `compiled_id`, falling back to interpreted R evaluation otherwise.
Compiled backends include formula-to-C++ compilation for ODE models, as
well as exact stochastic simulation algorithms (Gillespie direct,
Gibson-Bruck NRM, Anderson MNRM) for CTMC models.

For SSA solvers, `t_max` and `discard_transient` are in the same time
units as the propensity rates (typically dimensionless or in the natural
units of the model).

## See also

[`print.dyn_sim`](https://robustecologies.github.io/janos/reference/print.dyn_sim.md),
[`summary.dyn_sim`](https://robustecologies.github.io/janos/reference/summary.dyn_sim.md)
and
[`plot.dyn_sim`](https://robustecologies.github.io/janos/reference/plot.dyn_sim.md)
for the S3 methods operating on the returned object;
[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md)
for the constructor of the input model;
[`solver_rk45`](https://robustecologies.github.io/janos/reference/solver_rk45.md),
[`solver_euler_maruyama`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md)
and
[`solver_ssa_direct`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md)
as representative deterministic, stochastic-diffusive and
exact-stochastic solvers (see the solver family index for the full
list).

## Examples

``` r
if (FALSE) { # \dontrun{
# Lotka-Volterra predator-prey
model <- system_spec(
    rhs = list(
        prey ~ alpha * prey - beta * prey * predator,
        predator ~ delta * prey * predator - gamma * predator
    ),
    state_names = c("prey", "predator"),
    parms = list(alpha = 1.0, beta = 0.1, delta = 0.075, gamma = 1.5),
    init = c(prey = 40, predator = 9)
)
result <- dyn_sim(model, t_max = 50, solver = solver_rk4())
print(result)
summary(result)
plot(result, title = "Lotka-Volterra predator-prey")

# Simulate with adaptive RK45
result_rk45 <- dyn_sim(model, t_max = 50, solver = solver_rk45())
plot(result_rk45, title = "Lotka-Volterra (RK45)")

# SIR epidemic with Gillespie direct method
sir <- system_spec(
    stoichiometry = list(
        infection = c(S = -1L, I = 1L, R = 0L),
        recovery  = c(S = 0L,  I = -1L, R = 1L)
    ),
    propensities = list(
        infection ~ beta * S * I,
        recovery  ~ gamma * I
    ),
    state_names = c("S", "I", "R"),
    parms = list(beta = 0.001, gamma = 0.1),
    init = c(S = 999, I = 1, R = 0)
)
result_ssa <- dyn_sim(sir, t_max = 200, solver = solver_ssa_direct(),
                      discard_transient = 0)
plot(result_ssa, title = "SIR epidemic")
} # }
```
