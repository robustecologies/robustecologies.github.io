# Rosenbrock (Rodas3) implicit solver for stiff systems

Creates a solver specification for the Rodas3 linearly implicit
Rosenbrock method of order 4 with L-stability. The solver uses a
symbolically generated Jacobian (computed at compile time from the model
formulas) and adaptive step-size control via an embedded order-3 error
estimate. Suitable for stiff ODE systems where explicit methods require
impractically small time steps.

## Usage

``` r
solver_rosenbrock(
  dt_init = 0.01,
  atol = 1e-08,
  rtol = 1e-06,
  dt_min = 1e-10,
  dt_max = 1,
  max_steps = 1e+07,
  subsample = 0.1
)
```

## Arguments

- dt_init:

  Initial time step (default: 0.01).

- atol:

  Absolute error tolerance (default: 1e-8).

- rtol:

  Relative error tolerance (default: 1e-6).

- dt_min:

  Minimum allowed step size (default: 1e-10).

- dt_max:

  Maximum allowed step size (default: 1.0).

- max_steps:

  Maximum number of integration steps (default: 1e7).

- subsample:

  Approximate interval between stored output points (default: 0.1).

## Value

A list of class `solver_spec` with `method = "rosenbrock"`.

## Details

Rosenbrock methods solve a linear system (I/(gamma\*h) - J) \* k_i =
RHS_i at each stage, requiring one LU factorization per step but no
Newton iteration. The Jacobian J = df/dy is generated symbolically from
the model formulas, giving exact derivatives without finite-difference
approximation. The Rodas3 method uses 4 stages per step and is L-stable,
meaning it correctly damps very stiff components.

## References

Hairer, E. and Wanner, G. (1996). *Solving Ordinary Differential
Equations II: Stiff and Differential-Algebraic Problems*. Springer.
ISBN: 978-3-642-05221-7.

## See also

[`solver_rk45`](https://robustecologies.github.io/janos/reference/solver_rk45.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Van der Pol oscillator (stiff for large mu)
vdp <- system_spec(
    rhs = list(
        x ~ y,
        y ~ mu * (1 - x^2) * y - x
    ),
    state_names = c("x", "y"),
    parms = list(mu = 1000),
    init = c(x = 2, y = 0)
)
result <- dyn_sim(vdp, t_max = 10, solver = solver_rosenbrock())
print(result)
summary(result)
plot(result)
} # }
```
