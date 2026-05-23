# Fixed-step fourth-order Runge-Kutta solver

Creates a solver specification for the classical RK4 method with fixed
step size. This is the default solver and uses the same integration
scheme as the original janos implementation.

## Usage

``` r
solver_rk4(dt = 0.01, subsample = 10)
```

## Arguments

- dt:

  Integration time step (default: 0.01).

- subsample:

  Store every nth integration step to reduce memory (default: 10).

## Value

A list of class `solver_spec` with `method = "rk4"`.

## See also

[`solver_rk45`](https://robustecologies.github.io/janos/reference/solver_rk45.md),
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
result <- dyn_sim(model, t_max = 50, solver = solver_rk4())
} # }
```
