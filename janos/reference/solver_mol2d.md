# Specify a 2D method-of-lines solver

Creates a solver specification for integrating 2D PDE systems using the
method of lines with explicit RK4 time stepping. The spatial domain is
discretized on a Cartesian Nx x Ny grid and the resulting large ODE
system is advanced by a fixed-step fourth-order Runge-Kutta method.

## Usage

``` r
solver_mol2d(dt = 0.01, subsample = 10L)
```

## Arguments

- dt:

  Time step in simulation time units (default: 0.01).

- subsample:

  Store output every `subsample` steps (default: 10).

## Value

A `solver_spec` object with method `"mol2d"`.

## Details

The CFL stability condition for 2D diffusion with grid spacings dx, dy
and maximum diffusion coefficient D is dt \< dx^2 \* dy^2 / (2 \* D \*
(dx^2 + dy^2)). For dx = dy this simplifies to dt \< dx^2 / (4 \* D),
which is twice as restrictive as the 1D condition. A warning is issued
if the chosen dt exceeds this limit.

## See also

[`solver_mol`](https://robustecologies.github.io/janos/reference/solver_mol.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- model_spec(
    pde = list(u ~ D * lap(u)),
    state_names = "u",
    parms = list(D = 1.0),
    init = function(x, y) sin(pi * x) * sin(pi * y),
    spatial = list(
        domain = list(x = c(0, 1), y = c(0, 1)),
        n_grid = c(51, 51),
        bc = list(u = list(
            x = list(type = "dirichlet", left = 0, right = 0),
            y = list(type = "dirichlet", left = 0, right = 0)
        ))
    )
)
result <- dyn_sim(model, t_max = 0.1, solver = solver_mol2d(dt = 0.0001))
print(result)
plot(result)
} # }
```
