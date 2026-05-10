# Method of lines (MOL) solver for 1D PDE systems

Creates a solver specification for integrating 1D partial differential
equations via the method of lines (MOL). Spatial derivatives are
replaced with second-order central finite differences on a uniform grid,
converting the PDE into a system of coupled ODEs that is advanced in
time by the classical RK4 method.

## Usage

``` r
solver_mol(dt = 0.001, subsample = 10L)
```

## Arguments

- dt:

  Time step (default: 0.001).

- subsample:

  Store every nth integration step (default: 10).

## Value

A list of class `solver_spec` with `method = "mol"`.

## Details

The method of lines semi-discretizes the PDE by replacing spatial
derivatives with finite difference approximations on a grid of N points
with spacing dx = (x_max - x_min) / (N - 1). The resulting N coupled
ODEs per state variable are integrated with fixed-step RK4.

For pure diffusion (du/dt = D \* d2x(u)), the CFL stability condition
requires dt \< dx^2 / (2 \* D). A warning is issued if the
user-specified dt exceeds this bound.

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
heat <- model_spec(
    pde = list(u ~ D * d2x(u)),
    state_names = "u",
    parms = list(D = 0.01),
    spatial = list(
        domain = c(0, 1), n_grid = 101,
        bc = list(u = list(type = "dirichlet", left = 0, right = 0))
    ),
    init = function(x) sin(pi * x)
)
result <- dyn_sim(heat, t_max = 2, solver = solver_mol(dt = 0.001),
                  discard_transient = 0, verbose = FALSE)
print(result)
summary(result)
plot(result)
} # }
```
