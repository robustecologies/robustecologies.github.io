# Adaptive Dormand-Prince RK4(5) solver

Creates a solver specification for the Dormand-Prince embedded
Runge-Kutta method of order 4(5) with automatic step-size control. The
solver adapts the step size to maintain local truncation error within
the specified tolerances, eliminating the need for users to manually
select a fixed step.

## Usage

``` r
solver_rk45(
  dt_init = 0.01,
  atol = 1e-08,
  rtol = 1e-06,
  dt_min = 1e-06,
  dt_max = 1,
  max_steps = 1e+07,
  subsample = 0.1
)
```

## Arguments

- dt_init:

  Initial time step (default: 0.01). The solver adapts away from this
  value quickly.

- atol:

  Absolute error tolerance (default: 1e-8).

- rtol:

  Relative error tolerance (default: 1e-6).

- dt_min:

  Minimum allowed step size (default: 1e-6).

- dt_max:

  Maximum allowed step size (default: 1.0).

- max_steps:

  Maximum number of integration steps (default: 1e7). Prevents infinite
  loops if the solver cannot converge.

- subsample:

  Approximate interval between stored output points (default: 0.1). The
  solver stores output at evenly spaced times regardless of internal
  step sizes.

## Value

A list of class `solver_spec` with `method = "rk45"`.

## Details

The Dormand-Prince method uses six function evaluations per step (with
FSAL, "first same as last," optimization) to produce both a fourth-order
and a fifth-order solution. The difference between these estimates the
local truncation error, which is used to adjust the step size. The step
is accepted when the error satisfies:

\$\$\left\\\frac{y_5 - y_4}{atol + rtol \cdot \|y_4\|}\right\\\_\infty
\leq 1\$\$

When the error exceeds this bound, the step is rejected and retried with
a smaller step. After acceptance, the next step size is estimated from
the error ratio with a safety factor.

## References

Dormand, J. R., & Prince, P. J. (1980). A family of embedded Runge-Kutta
formulae. *Journal of Computational and Applied Mathematics*, 6(1),
19-26.
[doi:10.1016/0771-050X(80)90013-3](https://doi.org/10.1016/0771-050X%2880%2990013-3)

## See also

[`solver_rk4`](https://robustecologies.github.io/janos/reference/solver_rk4.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- system_spec(
    rhs = list(
        prey ~ alpha * prey - beta * prey * predator,
        predator ~ delta * prey * predator - gamma * predator
    ),
    state_names = c("prey", "predator"),
    parms = list(alpha = 1.0, beta = 0.1, delta = 0.075, gamma = 1.5),
    init = c(prey = 40, predator = 9)
)
result <- dyn_sim(model, t_max = 50, solver = solver_rk45())
print(result)
} # }
```
