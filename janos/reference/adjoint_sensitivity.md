# Adjoint sensitivity analysis for ODE models

Computes the gradient of a scalar objective function with respect to all
model parameters using the continuous adjoint method. The gradient is
obtained in a single backward pass, giving O(n_states + n_params) cost
regardless of the number of parameters.

## Usage

``` r
adjoint_sensitivity(
  model,
  objective = "terminal",
  objective_grad = NULL,
  t_max = 10,
  parms = NULL,
  init = NULL,
  n_checkpoints = 1000L,
  dt_backward = NULL,
  solver_forward = solver_rk45(),
  verbose = TRUE
)
```

## Arguments

- model:

  A `system_spec` with formula-based `rhs` and at least one parameter.

- objective:

  Character string or function specifying the objective. `"terminal"`
  (default) uses g = sum(y(T)) with dg/dy = 1 for all states. A custom
  function must accept a named numeric vector (terminal state) and
  return a scalar.

- objective_grad:

  Optional function returning the gradient dg/dy as a named numeric
  vector. Required when `objective` is a custom function. Ignored when
  `objective = "terminal"`.

- t_max:

  Simulation duration (default: 10).

- parms:

  Optional named list of parameter overrides.

- init:

  Optional named numeric vector of initial condition overrides.

- n_checkpoints:

  Number of evenly spaced checkpoints to store from the forward
  trajectory for interpolation during backward integration (default:
  1000).

- dt_backward:

  Step size for backward RK4 integration. If NULL (default), set to half
  the forward trajectory time spacing.

- solver_forward:

  Solver for the forward pass (default:
  [`solver_rk45()`](https://robustecologies.github.io/janos/reference/solver_rk45.md)).

- verbose:

  Logical; print progress messages (default: TRUE).

## Value

An S3 object of class `sensitivity_result` containing:

- gradient:

  Named numeric vector of dg/dp for each parameter

- relative_sensitivity:

  Named numeric vector of elasticities: (dg/dp) \* p / g, measuring the
  proportional response of the objective to a proportional change in
  each parameter

- objective_value:

  Scalar value of g(y(T))

- forward_trajectory:

  Data frame of the forward trajectory

- lambda_trajectory:

  Data frame of the adjoint variable trajectory

- model:

  The system_spec used

- parms:

  Named list of parameter values

- t_max:

  Simulation duration

## Details

For an ODE system dy/dt = f(y, p) with scalar objective g(y(T)), the
adjoint method proceeds in three steps. First, the forward ODE is
integrated to the terminal time T, storing the trajectory at
checkpoints. Second, the adjoint variable lambda(t) is integrated
backward from the terminal condition lambda(T) = dg/dy evaluated at t=T
via the adjoint ODE d(lambda)/dt = -J(y(t))^T \* lambda(t), where J is
the Jacobian of f with respect to y. Third, the parameter gradient is
accumulated by quadrature: dg/dp = -integral from 0 to T of (df/dp)^T \*
lambda(t) dt.

Both the transposed Jacobian and the df/dp matrix are generated
symbolically from the model formulas and compiled to C++ for
high-performance evaluation. The forward trajectory is interpolated
linearly between checkpoints during backward integration.

## References

Pontryagin, L. S., Boltyanskii, V. G., Gamkrelidze, R. V., & Mishchenko,
E. F. (1962). *The Mathematical Theory of Optimal Processes*.
Wiley-Interscience.

Cao, Y., Li, S., Petzold, L., & Serban, R. (2003). Adjoint sensitivity
analysis for differential-algebraic equations: the adjoint DAE system
and its numerical solution. *SIAM Journal on Scientific Computing*,
24(3), 1076-1089.
[doi:10.1137/S1064827501380630](https://doi.org/10.1137/S1064827501380630)

## See also

[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
[`equilibrium`](https://robustecologies.github.io/janos/reference/equilibrium.md),
[`continuation`](https://robustecologies.github.io/janos/reference/continuation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Exponential decay: dy/dt = -k * y
decay <- system_spec(
    rhs = list(y ~ -k * y),
    state_names = "y",
    parms = list(k = 0.5),
    init = c(y = 1.0)
)

# Sensitivity of terminal state to k
sens <- adjoint_sensitivity(decay, t_max = 2)
print(sens)
summary(sens)
plot(sens)
plot(sens, type = "trajectory")

# Compare with analytical: dg/dk = -y0 * T * exp(-k * T)
analytical <- -1.0 * 2.0 * exp(-0.5 * 2.0)
print(paste("Adjoint:", sens$gradient["k"], "Analytical:", analytical))
} # }
```
